# glue_F2_over_F_vs_gsq_claude.py
# F^2 (0++) / F (l=1) glueball mass ratio vs g^2, L=1,2.
#   R = Delta_{F^2,0++} / Delta_{F,l=1}
# Both from the GEVP: F^2 0++ = fit_perm SINGLE state "0" over t=[0.2,0.4] (prefix gevp_f2);
#   F l=1 = fit_perm variance-avg over t=[0.2,1.0] (prefix gevp_F).  Naive error propagation.
# Free single-particle refs on S2: F^2 0++ = 2 sqrt2 (2 photons), F l=1 = sqrt2 -> R_free = 2.
import subprocess, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NFS = [2, 4, 6]
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0]}
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
Lmk = {1: "o", 2: "s"}
FREE = 2.0
F2WIN = (0.2, 0.4)
FWIN = (0.2, 1.0)


def glue_mass(pref, tag, states, tlo, thi):
    jk = "%s_%s_jk_claude.dat" % (pref, tag)
    try:
        open(jk).close()
    except OSError:
        return None, None
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(tlo), str(thi), states],
                         capture_output=True, text=True).stdout
    key = "variance-avg" if "," in states else ("state %s:" % states)
    for line in out.splitlines():
        if key in line:
            seg = line.split("M =")[-1] if "==>" in line else line.split("M=")[-1]
            m = seg.strip().split()[0]
            try:
                v = float(m.split("(")[0])
                e = float(m.split("(")[1].rstrip(")"))
            except (ValueError, IndexError):
                return None, None
            return (v, e) if np.isfinite(v) else (None, None)
    return None, None


rows = []   # (L, nf, g, r, er)
for L in [1, 2]:
    for nf in NFS:
        for g in GS[L]:
            if abs(g - 3.0) < 1e-9:      # omit g^2=3 (L2 g3.0: 0++ too noisy, NM 2026-08-18)
                continue
            tag = "Nf%d_g%.6f_L%d" % (nf, g, L)
            m2, e2 = glue_mass("gevp_f2", tag, "0", F2WIN[0], F2WIN[1])
            mF, eF = glue_mass("gevp_F", tag, "0,1,2", FWIN[0], FWIN[1])
            if m2 is None or mF is None:
                continue
            r = m2 / mF
            er = r * math.sqrt((e2 / m2) ** 2 + (eF / mF) ** 2)
            rows.append((L, nf, g, r, er))
            print("  L%d Nf%d g%.1f : F2/F = %.4f(%.4f)" % (L, nf, g, r, er))

fig, ax = plt.subplots(figsize=(7.5, 5.5))
for L in [1, 2]:
    for nf in NFS:
        d = sorted([x for x in rows if x[0] == L and x[1] == nf], key=lambda x: x[2])
        if not d:
            continue
        ax.errorbar([x[2] for x in d], [x[3] for x in d], [x[4] for x in d],
                    marker=Lmk[L], color=nfcol[nf], capsize=3, ls="-", lw=0.8,
                    mfc="none" if L == 2 else nfcol[nf], label="Nf%d L%d" % (nf, L))
# free reference removed (NM 2026-08-18: 2 is not the correct free value here)
ax.set_xlabel(r"$g_0^2$")
ax.set_ylabel(r"$\Delta_{F^2,0^{++}}/\Delta_{F,\ell=1}$")
ax.set_title(r"$F^2$ ($0^{++}$) / $F$ ($\ell{=}1$) vs $g_0^2$  (L1 filled, L2 open)")
ax.grid(alpha=0.3)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig("figs/ratio_F2_over_F_vs_gsq_claude.png", dpi=150)
print("# wrote ratio_F2_over_F_vs_gsq_claude.png")

out = ["# F^2 (0++) / F(l=1) mass ratio vs g^2, L=1,2 (both GEVP).  Free ref 2.", "",
       "| L | gsq | Nf | F2_0++/F1 |", "|---|-----|----|-----------|"]
for L, nf, g, r, er in sorted(rows):
    out.append("| %d | %.1f | %d | %.4f(%.4f) |" % (L, g, nf, r, er))
open("ratio_F2_over_F_vs_gsq_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote ratio_F2_over_F_vs_gsq_claude.md")
