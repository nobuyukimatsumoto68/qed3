# ratio_axial_tp_l2l1_vs_a2_claude.py
# The DIMENSIONLESS axial tp ratio  m0(ell=2)/m0(ell=1)  vs lattice spacing squared a^2 ~ 1/L^2.
# READS the consolidated mass summary  axial_tp_masses_summary_claude.md  (written by
# effmass_axial_tp_expfit_claude.py) -- single source of truth, canonical const+exp tp windows.
# Ratio r = m0(l2)/m0(l1) per (L, gsq, Nf); the a_t normalization cancels in the ratio.
# NOTE: the summary carries only (mass, error), not the per-sample covariance, so the error here is
# NAIVE (uncorrelated) propagation  dr/r = sqrt((e2/m2)^2 + (e1/m1)^2).  (l2 and l1 are positively
# correlated across configs, so this is a conservative upper bound on the true ratio error.)
# x-axis a^2 ~ 1/L^2 (a/R ~ 1/L on the refined sphere).  Free/CFT axial ratio Delta_2/Delta_1 = 3/2.
import re, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SUMMARY = "axial_tp_masses_summary_claude.md"
NCFG_MIN = 100   # global cut: omit any ensemble point with fewer than 100 configs (NM 2026-08-18)
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
FREE = 1.5   # free/CFT axial ratio Delta_(l=2)/Delta_(l=1) = 3/2  (Delta_l = l+1)
dxNf = {2: -0.012, 4: 0.0, 6: 0.012}


# ---------- parse the summary table ----------
# cols: ell | irrep | L | 1/L^2 | gsq | Nf | fit dt lo | fit dt hi | a_t m0 | err | ncfg
M = {}      # (ell, L, g, nf) -> (m0, err, ncfg)
invL2 = {}  # L -> 1/L^2
for line in open(SUMMARY):
    s = line.strip()
    if not s.startswith("|"):
        continue
    c = [x.strip() for x in s.strip("|").split("|")]
    if len(c) != 11 or c[0] in ("ell", "-----") or c[0].startswith("-"):
        continue
    try:
        ell = int(c[0])
        L = int(c[2])
        iL2 = float(c[3])
        g = float(c[4])
        nf = int(c[5])
        m0 = float(c[8])
        err = float(c[9])
        ncfg = int(c[10])
    except ValueError:
        continue
    M[(ell, L, g, nf)] = (m0, err, ncfg)
    invL2[L] = iL2

# ---------- form l2/l1 ratios ----------
ratios = []   # (L, invL2, g, nf, r, er, ncfg)
for (ell, L, g, nf), (m2, e2, n2) in M.items():
    if ell != 2:
        continue
    if (1, L, g, nf) not in M:
        continue
    m1, e1, n1 = M[(1, L, g, nf)]
    if m1 <= 0 or m2 <= 0:
        continue
    if min(n1, n2) < NCFG_MIN:     # global ncfg cut
        continue
    r = m2 / m1
    er = r * math.sqrt((e2 / m2) ** 2 + (e1 / m1) ** 2)
    ratios.append((L, invL2[L], g, nf, r, er, min(n1, n2)))
    print("  L%d g%.1f Nf%d : ratio l2/l1 = %.4f(%.4f)  ncfg=%d" % (L, g, nf, r, er, min(n1, n2)))

# ---------- plot vs 1/L^2 ----------
fig, ax = plt.subplots(figsize=(7.5, 5.5))
for nf in NFS:
    first = True
    for L, iL2, g, nf_, r, er, n in sorted(ratios):
        if nf_ != nf:
            continue
        ax.errorbar(iL2 + dxNf[nf], r, yerr=er, marker="o", ms=6, capsize=3,
                    color=nfcol[nf], ls="none", label=("Nf%d" % nf) if first else None)
        first = False
ax.axhline(FREE, color="gray", ls="--", lw=1.0)
ax.text(0.02, FREE, r" free $3/2$", color="gray", fontsize=9, va="bottom",
        transform=ax.get_yaxis_transform())
ax.set_xlabel(r"$a^2 \sim 1/L^2$")
ax.set_ylabel(r"$m_0(\ell{=}2)/m_0(\ell{=}1)$  (axial tp)")
ax.set_title(r"axial tp $\ell{=}2/\ell{=}1$ vs $a^2$ (const+exp $m_0$ from summary)")
ax.set_xlim(-0.05, 1.1)
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("figs/ratio_axial_tp_l2l1_vs_a2_claude.png", dpi=150)
print("# wrote ratio_axial_tp_l2l1_vs_a2_claude.png")

# ---------- md dump ----------
out = ["# Axial tp ratio m0(l=2)/m0(l=1) vs a^2~1/L^2 (from axial_tp_masses_summary_claude.md)", "",
       "Naive (uncorrelated) error propagation from the summary masses/errors.", "",
       "| L | 1/L^2 | gsq | Nf | ratio l2/l1 | ncfg |",
       "|---|-------|-----|----|-------------|------|"]
for L, iL2, g, nf, r, er, n in sorted(ratios):
    out.append("| %d | %.4f | %.1f | %d | %.4f(%.4f) | %d |" % (L, iL2, g, nf, r, er, n))
open("ratio_axial_tp_l2l1_vs_a2_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote ratio_axial_tp_l2l1_vs_a2_claude.md")
