# glue_F2_over_axial_vs_gsq_claude.py
# F^2 (0++) glueball / axial(l=1) mass ratio vs g^2, L=1,2.
#   R = Delta_{F^2,0++} / Delta_{A,l=1}
# F^2 0++ = GEVP (F^2-v2 basis, F^4-included) fit_perm SINGLE state "0" over t=[0.2,0.4]
#   (prefix gevp_f2, PHYSICAL units).  axial l=1 = const+exp m0 from axial_tp_masses_summary_claude.md.
# Free single-particle refs on S2: F^2 0++ = 2 sqrt2 (2 photons), axial l=1 = 2 -> R_free = sqrt2 ~ 1.414.
import subprocess, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NFS = [2, 4, 6]
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0]}     # L=1,2 first
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
Lmk = {1: "o", 2: "s"}
SUMMARY = "axial_tp_masses_summary_claude.md"
FREE = math.sqrt(2.0)
F2WIN = (0.2, 0.4)


def axial_summary():
    d = {}
    for line in open(SUMMARY):
        s = line.strip()
        if not s.startswith("|"):
            continue
        c = [x.strip() for x in s.strip("|").split("|")]
        if len(c) != 11 or not c[0].isdigit() or int(c[0]) != 1:
            continue
        d[(int(c[2]), float(c[4]), int(c[5]))] = (float(c[8]), float(c[9]))
    return d


axsum = axial_summary()


def f2_mass(tag):
    # single-state "0" fit_perm on the 0++ jk dump
    jk = "gevp_f2_%s_jk_claude.dat" % tag
    try:
        open(jk).close()
    except OSError:
        return None, None
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(F2WIN[0]), str(F2WIN[1]), "0"],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "state 0:" in line:
            seg = line.split("M=")[-1]
            m = seg.strip().split()[0]
            try:
                v = float(m.split("(")[0])
                e = float(m.split("(")[1].rstrip(")"))
            except (ValueError, IndexError):
                return None, None
            if not np.isfinite(v):
                return None, None
            return v, e
    return None, None


rows = []   # (L, nf, g, r, er)
for L in [1, 2]:
    for nf in NFS:
        for g in GS[L]:
            if (L, g, nf) not in axsum:
                continue
            tag = "Nf%d_g%.6f_L%d" % (nf, g, L)
            m2, e2 = f2_mass(tag)
            if m2 is None:
                continue
            a, ea = axsum[(L, g, nf)]
            r = m2 / a
            er = r * math.sqrt((e2 / m2) ** 2 + (ea / a) ** 2)
            rows.append((L, nf, g, r, er))
            print("  L%d Nf%d g%.1f : F2/A1 = %.4f(%.4f)" % (L, nf, g, r, er))

fig, ax = plt.subplots(figsize=(7.5, 5.5))
for L in [1, 2]:
    for nf in NFS:
        d = sorted([x for x in rows if x[0] == L and x[1] == nf], key=lambda x: x[2])
        if not d:
            continue
        ax.errorbar([x[2] for x in d], [x[3] for x in d], [x[4] for x in d],
                    marker=Lmk[L], color=nfcol[nf], capsize=3, ls="-", lw=0.8,
                    mfc="none" if L == 2 else nfcol[nf], label="Nf%d L%d" % (nf, L))
# free reference removed (NM 2026-08-18: sqrt2 is not the correct free value here)
ax.set_xlabel(r"$g_0^2$")
ax.set_ylabel(r"$\Delta_{F^2,0^{++}}/\Delta_{A,\ell=1}$")
ax.set_title(r"$F^2$ ($0^{++}$) / axial $\ell{=}1$ vs $g_0^2$  (L1 filled, L2 open)")
ax.grid(alpha=0.3)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig("figs/ratio_F2_over_axial_vs_gsq_claude.png", dpi=150)
print("# wrote ratio_F2_over_axial_vs_gsq_claude.png")

out = ["# F^2 (0++) / axial(l=1) mass ratio vs g^2 (axial from const+exp summary), L=1,2", "",
       "0++ = gevp_f2 fit_perm state 0 over t[0.2,0.4].  Free ref sqrt2 ~ 1.414.", "",
       "| L | gsq | Nf | F2_0++/A1 |", "|---|-----|----|-----------|"]
for L, nf, g, r, er in sorted(rows):
    out.append("| %d | %.1f | %d | %.4f(%.4f) |" % (L, g, nf, r, er))
open("ratio_F2_over_axial_vs_gsq_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote ratio_F2_over_axial_vs_gsq_claude.md")
