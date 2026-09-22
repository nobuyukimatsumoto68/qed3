#!/usr/bin/env python3
# t00_stress_gevp_claude.py  [CHUNK 1b: multi-displacement variational GEVP for fermionic T_00]
#
# Extends the single-operator T_00 two-point (t00_stress_claude.py) to a variational basis of temporal
# displacements Dt = 1,2,3.  The single Dt=1 operator only reaches the free target 3/R=0.567 near the
# window edge (effmass approaches from above); a GEVP over {Dt} projects the \Delta=3 ground state earlier.
#
# Run (free-limit, single config, exact):
#   ENS=free NVDIR=distill_Nv24 LREF=1 DISPS=1,2,3 T0=3 python3 t00_stress_gevp_claude.py
#
# Operator family (ell=0, 0^{++}, all in the T_00 channel):
#   O_{Dt}(t) = sum_x w_x psibar(x,t) sigma_3 (1/2)[ psi(x,t+Dt) - psi(x,t-Dt) ] ,  w_x = A_x Y00 .
# Cross two-point (sink displacement Da at s+dt, source displacement Db at s; translation-avg over s):
#   C_{ab}(dt) = -(1/4) sum_{sx,sy=+-1} sx sy sum_{x,y} w_x w_y
#                  Tr[ sigma_3 S(x,s+dt+sx*Da; y,s) sigma_3 S(y,s+sy*Db; x,s+dt) ] .
# S = AblkS (fs_gevp_point make_config); contact already subtracted in the perambulator (A - 1/2 I at eq-time).
# See t00_stress_impl_plan_claude.md.

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
DTMAX = int(os.environ.get("DTMAX", "24"))
T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "1"))
NCFG = int(os.environ.get("NCFG", "0"))
DISPS = [int(x) for x in os.environ.get("DISPS", "1,2,3").split(",")]
ROVER = 1.0 / 0.189                                # free 1/R (a_t m): sigma=2/R, T00=3/R, 2m=4/R

SIG3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


def matrix_one_config(k, w, disps):
    # Connected T_00 GEVP matrix C[a,b,dt] for one config (a=sink disp, b=source disp).
    AblkS, AblkSt, twin, nsite, U, tau = fg.make_config(k)
    nop = len(disps)
    C = np.full((nop, nop, DTMAX), np.nan)
    signs = (1, -1)
    for dt in range(DTMAX):
        for a in range(nop):
            Da = disps[a]
            for b in range(nop):
                Db = disps[b]
                # all four block time-args in [0, twin-1]
                s_lo = max(Db, Da - dt)
                s_hi = min(twin - 1 - Db, twin - 1 - dt - Da)
                if s_hi < s_lo:
                    continue
                acc = 0.0
                cnt = 0
                for s in range(s_lo, s_hi + 1):
                    val = 0.0
                    for sx in signs:
                        B1 = AblkS(s + dt + sx * Da, s)          # [x,a,y,b]
                        for sy in signs:
                            B2 = AblkS(s + sy * Db, s + dt)      # [y,c,x,d]
                            t = np.einsum("ab,pbqc,cd,qdpa,p,q->", SIG3, B1, SIG3, B2, w, w, optimize=True)
                            val += (sx * sy) * t
                    acc += val
                    cnt += 1
                C[a, b, dt] = (-0.25 * acc / cnt).real
    return C


def main():
    tag = dc.ENS
    dual = dc.dual_areas_from_mesh()
    w = dual * dc.Y00
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    if not ks:
        print("# no perambulators found in %s" % dc.PERAM_DIR)
        return
    nop = len(DISPS)
    print("# ENS=%s NVDIR=%s  nsite=%d ncfg=%d  DISPS=%s  T0=%d DTMAX=%d"
          % (tag, dc.NVDIR, dual.shape[0], len(ks), DISPS, T0, DTMAX))

    allC = np.array([matrix_one_config(k, w, DISPS) for k in ks])   # (ncfg, nop, nop, DTMAX)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))                   # symmetrize
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])

    def gevp_effmass(Cm):
        ev = np.full((DTMAX, nop), np.nan)
        for dt in range(DTMAX):
            if np.any(~np.isfinite(Cm[:, :, dt])) or np.any(~np.isfinite(Cm[:, :, T0])):
                continue
            try:
                ev[dt] = fg.gevp(Cm[:, :, dt], Cm[:, :, T0])
            except Exception:
                pass
        with np.errstate(all="ignore"):
            return np.log(ev[:-1] / ev[1:])

    em_c = gevp_effmass(blk.mean(0))
    if nb < 2:
        em_err = np.zeros_like(em_c)                 # single config (free limit): exact
    else:
        ems = np.array([gevp_effmass(np.delete(blk, i, 0).mean(0)) for i in range(nb)])
        em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))

    print("\n#  dt |  " + "  ".join("m%d(err)   " % n for n in range(nop))
          + "  [free refs sigma=%.3f T00=%.3f 2m=%.3f]" % (2.0 / ROVER, 3.0 / ROVER, 4.0 / ROVER))
    for dt in range(T0, min(em_c.shape[0], 22)):
        row = "  ".join("%7.4f(%.4f)" % (em_c[dt, n], em_err[dt, n]) for n in range(nop))
        print("#  %2d | %s" % (dt, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(em_c.shape[0])
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    for val, lab, col in ((2.0 / ROVER, r"$\sigma=2/R$", "tab:green"),
                          (3.0 / ROVER, r"$T_{00}=3/R$", "tab:red"),
                          (4.0 / ROVER, r"$2m=4/R$", "gray")):
        ax.axhline(val, color=col, ls="--", lw=1, alpha=0.7)
        ax.text(T0 + 0.1, val + 0.006, lab, fontsize=9, color=col)
    cols = ["tab:red", "tab:blue", "tab:purple"]
    mkr = ["o", "s", "^"]
    for n in range(nop):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n % 3], marker=mkr[n % 3], ms=5,
                    lw=1.1, capsize=2.5, label="GEVP state %d" % n)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(em_c.shape[0], 22))
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"$T_{00}$ GEVP  DISPS=%s  T0=%d  %s  %d cfg" % (DISPS, T0, tag, ncfg), fontsize=10.5)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_stress_gevp_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
