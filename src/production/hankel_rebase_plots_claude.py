#!/usr/bin/env python3
# hankel_rebase_plots_claude.py
#   Produce labeled effmass figures for the shortlisted (and rejected) block-Hankel + rebase
#   configurations, using the cached store.  Figures saved to figs/ with clear names.
#   2m_{PS} ~ 0.64 (m_{PS} ~ 0.32) drawn as a reference line for the two-meson level.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import hankel_rebase_scan_claude as hs

store, twin, ncfg, tag = hs.build_store()
M2PS = 0.644     # 2 m_{PS}, m_{PS} ~ 0.322
os.makedirs("figs", exist_ok=True)


def make_fig(fname, title, nooa, offsets, stages, t0, note):
    em_c, em_err, ems, nbin = hs.eval_config(store, twin, nooa, offsets, stages, t0)
    tmax = em_c.shape[0]
    nk = em_c.shape[1]
    ts = np.arange(tmax)
    cols = ["tab:green", "tab:red", "tab:blue"]
    mkr = ["o", "s", "^"]
    labs = ["state 0 (one-meson-rich ground)", "state 1", "state 2 (two-meson)"]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.axhline(M2PS, color="gray", ls="--", lw=1, alpha=0.7)
    ax.text(tmax * 0.62, M2PS + 0.012, r"$2m_{PS}\approx0.644$", fontsize=9, color="gray")
    for n in range(nk):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.2)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n % 3], marker=mkr[n % 3],
                    ms=5, lw=1.1, capsize=2.5, label=labs[n] if n < 3 else "state %d" % n)
    for (tr, nk_s) in stages:
        ax.axvline(tr, color="k", ls=":", lw=0.9, alpha=0.5)
    ax.set_ylim(0.2, 0.9)
    ax.set_xlim(t0, min(tmax, 22))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(title, fontsize=11)
    ax.text(0.02, 0.02, note, transform=ax.transAxes, fontsize=8, va="bottom",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = "figs/" + fname
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("# -> %s" % out)


# ---- shortlist ----
make_fig("hankel_scan_pick1_Dt12_2op_T0-3_reb5_claude.png",
         "PICK 1  Dt=1,2 [0,1,2]  2-op  T0=3  rebase 2@t=5  (baseline sweet spot)",
         1, [0, 1, 2], [(5, 2)], 3,
         "m0=0.470(10)  m1=0.644(3) FLAT (slope -0.003(0))\ntwo-meson right at 2mPS; best all-around, range to t~20")

make_fig("hankel_scan_pick2_Dt3_2op_T0-3_reb5_claude.png",
         "PICK 2  Dt=3 [0,3]  2-op  T0=3  rebase 2@t=5  (cheap 2-block alt)",
         1, [0, 3], [(5, 2)], 3,
         "m0=0.479(9)  m1=0.648(3) FLAT (slope -0.004(0))\n2 blocks only; slightly tighter m0 at mid-t")

make_fig("hankel_scan_pick3_Dt12_2op_T0-2_reb5_claude.png",
         "PICK 3  Dt=1,2 [0,1,2]  2-op  T0=2  rebase 2@t=5  (earliest onset)",
         1, [0, 1, 2], [(5, 2)], 2,
         "m0=0.483(9)  m1=0.649(3) FLAT (slope -0.004(0))\nT0=2: earliest plateau, tightest m1 S/N")

make_fig("hankel_scan_pick4_Dt24_2op_T0-3_reb5_claude.png",
         "PICK 4  Dt=2,4 [0,2,4]  2-op  T0=3  rebase 2@t=5  (early m1 plateau)",
         1, [0, 2, 4], [(5, 2)], 3,
         "m0=0.430(12)  m1=0.633(3) FLAT (slope -0.002(1))\nfastest m1 plateau (t~4) but m1 slightly under 2mPS, range to t~13")

make_fig("hankel_scan_pick5_Dt12_3op_T0-3_nkeep3_claude.png",
         "PICK 5  Dt=1,2 [0,1,2]  3-op(+O_A)  T0=3  rebase 3@t=5  (resolves 3 levels)",
         0, [0, 1, 2], [(5, 3)], 3,
         "m0=0.438(22)  m1~0.48 (O_A state, tilts)  m2=0.637(3) FLAT two-meson\nuse if the intermediate O_A state is wanted; two-meson = level 2")

# ---- rejects ----
make_fig("hankel_scan_REJECT_Dt1_2op_tilt_claude.png",
         "REJECT  Dt=1 [0,1]  2-op  T0=3  rebase 2@t=5  (m1 TILTS)",
         1, [0, 1], [(5, 2)], 3,
         "m1 slope -0.0097(3): NOT flat -> under-cleaned, residual contamination in state 1")

make_fig("hankel_scan_REJECT_Dt48_2op_overclean_claude.png",
         "REJECT  Dt=4,8 [0,4,8]  2-op  T0=3  rebase 2@t=5  (over-cleaned)",
         1, [0, 4, 8], [(5, 2)], 3,
         "m1 flat but pulled to 0.612 (< 2mPS); m0 has NO plateau (dips then rises); range only t<=13")

print("# all figures written to figs/")
