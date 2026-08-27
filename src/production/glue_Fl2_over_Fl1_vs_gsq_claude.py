# glue_Fl2_over_Fl1_vs_gsq_claude.py
# F(l=2, H) / F(l=1) glueball mass ratio vs g^2, L=1,2,3 (both from GEVP fit_perm variance-avg).
#   R = Delta_{F,l=2} / Delta_{F,l=1}
# F l=1 = fit_perm [0.2,1.0] (gevp_F); F l=2 = fit_perm [0.2,0.6] (gevp_Fl2).  Naive error propagation.
# Free single-particle refs on S2: F l=2 = sqrt3, F l=1 = sqrt2 -> R_free = sqrt(3/2) ~ 1.225.
import subprocess, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NFS = [2, 4, 6]
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5]}
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
Lmk = {1: "o", 2: "s", 3: "D"}
FREE = math.sqrt(1.5)
F1WIN = (0.2, 1.0)
F2WIN = (0.2, 0.6)


def glue_mass(pref, tag, tlo, thi):
    jk = "%s_%s_jk_claude.dat" % (pref, tag)
    try:
        open(jk).close()
    except OSError:
        return None, None
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(tlo), str(thi), "0,1,2"],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "variance-avg" in line:
            seg = line.split("M =")[-1] if "==>" in line else line.split("M=")[-1]
            m = seg.strip().split()[0]
            try:
                v = float(m.split("(")[0])
                e = float(m.split("(")[1].rstrip(")"))
            except (ValueError, IndexError):
                return None, None
            return (v, e) if np.isfinite(v) else (None, None)
    return None, None


rows = []
for L in [1, 2, 3]:
    for nf in NFS:
        for g in GS[L]:
            tag = "Nf%d_g%.6f_L%d" % (nf, g, L)
            m1, e1 = glue_mass("gevp_F", tag, F1WIN[0], F1WIN[1])
            m2, e2 = glue_mass("gevp_Fl2", tag, F2WIN[0], F2WIN[1])
            if m1 is None or m2 is None:
                continue
            r = m2 / m1
            er = r * math.sqrt((e2 / m2) ** 2 + (e1 / m1) ** 2)
            rows.append((L, nf, g, r, er))
            print("  L%d Nf%d g%.1f : Fl2/Fl1 = %.4f(%.4f)" % (L, nf, g, r, er))

fig, ax = plt.subplots(figsize=(7.5, 5.5))
for L in [1, 2, 3]:
    for nf in NFS:
        d = sorted([x for x in rows if x[0] == L and x[1] == nf], key=lambda x: x[2])
        if not d:
            continue
        ax.errorbar([x[2] for x in d], [x[3] for x in d], [x[4] for x in d],
                    marker=Lmk[L], color=nfcol[nf], capsize=3, ls="-", lw=0.8,
                    mfc="none" if L == 2 else nfcol[nf], label="Nf%d L%d" % (nf, L))
ax.axhline(FREE, color="gray", ls="--", lw=1.0)
ax.text(0.02, FREE, r" free $\sqrt{3/2}$", color="gray", fontsize=9, va="bottom",
        transform=ax.get_yaxis_transform())
ax.set_xlabel(r"$g_0^2$")
ax.set_ylabel(r"$\Delta_{F,\ell=2}/\Delta_{F,\ell=1}$")
ax.set_title(r"$F$ ($\ell{=}2$, H) / $F$ ($\ell{=}1$) vs $g_0^2$  (L1 o, L2 open, L3 D)")
ax.grid(alpha=0.3)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig("figs/ratio_Fl2_over_Fl1_vs_gsq_claude.png", dpi=150)
print("# wrote ratio_Fl2_over_Fl1_vs_gsq_claude.png")

out = ["# F(l=2)/F(l=1) mass ratio vs g^2, L=1,2,3 (both GEVP).  Free ref sqrt(3/2) ~ 1.225.", "",
       "| L | gsq | Nf | Fl2/Fl1 |", "|---|-----|----|---------|"]
for L, nf, g, r, er in sorted(rows):
    out.append("| %d | %.1f | %d | %.4f(%.4f) |" % (L, g, nf, r, er))
open("ratio_Fl2_over_Fl1_vs_gsq_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote ratio_Fl2_over_Fl1_vs_gsq_claude.md")
