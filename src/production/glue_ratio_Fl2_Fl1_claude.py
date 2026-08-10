# glue_ratio_Fl2_Fl1_claude.py
# Ratio of glueball F effective masses Delta_{F,l=2}/Delta_{F,l=1} for the massless redo
# ensembles (per-m GEVP, kmin=20).  F l=1 = T1 (3 m-grounds), fit t in [0.2,1.0];
# F l=2 = H (5 m-grounds), fit t in [0.2,0.6].  Masses from fit_perm variance-avg over the
# per-m grounds.  Free single-photon reference: omega_l = sqrt(l+1) (l=1 -> sqrt2, l=2 -> sqrt3) -> ratio sqrt(3/2).
# Two figures: vs gsq (L1,L2 markers) and vs 1/L^2 (per-Nf x-displacement).
import subprocess, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NFS = [2, 4, 6]
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
Lmk = {1: "o", 2: "s", 3: "D", 4: "^"}
# channel -> (dat/jk prefix, states list, fit window)
CH = {
    "Fl1": ("gevp_F",   "0,1,2",     (0.2, 1.0)),
    "Fl2": ("gevp_Fl2", "0,1,2,3,4", (0.2, 0.6)),
}


def fit(jk, states, tlo, thi):
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(tlo), str(thi), states],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "variance-avg" in line:
            seg = line.split("M =")[-1]
            m = seg.strip().split()[0]
            return float(m.split("(")[0]), float(m.split("(")[1].rstrip(")"))
    return None, None


rows = []
for L in [1, 2, 3, 4]:
    for nf in NFS:
        for g in GS[L]:
            tag = "Nf%d_g%.6f_L%d" % (nf, g, L)
            vals = {}
            for ch, (pref, states, (tlo, thi)) in CH.items():
                jk = "%s_%s_jk_claude.dat" % (pref, tag)
                M, e = fit(jk, states, tlo, thi)
                vals[ch] = (M, e)
            (m1, e1), (m2, e2) = vals["Fl1"], vals["Fl2"]
            if m1 is None or m2 is None:
                continue
            r = m2 / m1
            er = r * math.sqrt((e1 / m1) ** 2 + (e2 / m2) ** 2)
            rows.append((L, nf, g, m1, e1, m2, e2, r, er))

# ---- table ----
print("# L Nf gsq  m_Fl1        m_Fl2        ratio")
for L, nf, g, m1, e1, m2, e2, r, er in rows:
    print("%d %d %.1f  %.4f(%.4f) %.4f(%.4f)  %.4f(%.4f)" % (L, nf, g, m1, e1, m2, e2, r, er))

# ---- vs gsq ----
fig, ax = plt.subplots(figsize=(7, 5))
for L in [1, 2, 3, 4]:
    for nf in NFS:
        d = [x for x in rows if x[0] == L and x[1] == nf]
        if not d:
            continue
        gg = [x[2] for x in d]
        rr = [x[7] for x in d]
        ee = [x[8] for x in d]
        ax.errorbar(gg, rr, ee, marker=Lmk[L], color=nfcol[nf], capsize=3, ls="-", lw=0.8,
                    mfc="none" if L == 2 else nfcol[nf],
                    label="Nf%d L%d" % (nf, L))
ax.axhline(math.sqrt(1.5), color="gray", ls="--", lw=1.0)
ax.text(0.05, math.sqrt(1.5) + 0.02, r"free $\sqrt{3/2}$", color="gray", fontsize=9,
        transform=ax.get_yaxis_transform())
ax.set_xlabel(r"$g_0^2$")
ax.set_ylabel(r"$\Delta_{F,\ell=2}/\Delta_{F,\ell=1}$")
ax.set_title(r"Glueball $F$ ratio $\ell=2/\ell=1$ vs $g_0^2$")
ax.grid(alpha=0.3)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig("figs/glue_ratio_Fl2_Fl1_vs_gsq_claude.png", dpi=150)
print("wrote glue_ratio_Fl2_Fl1_vs_gsq_claude.png")

# ---- vs 1/L^2 (per-Nf x-displacement) ----
fig, ax = plt.subplots(figsize=(7, 5))
invL2 = {1: 1.0, 2: 0.25, 3: 1.0/9.0, 4: 0.0625}
dx = {2: -0.012, 4: 0.0, 6: 0.012}
for nf in NFS:
    for g in GS[1] + GS[2] + GS[3] + GS[4]:
        d = [x for x in rows if x[1] == nf and abs(x[2] - g) < 1e-9]
        if not d:
            continue
        xx = [invL2[x[0]] + dx[nf] for x in d]
        rr = [x[7] for x in d]
        ee = [x[8] for x in d]
        ax.errorbar(xx, rr, ee, marker="o", color=nfcol[nf], capsize=3, ls="none",
                    label="Nf%d" % nf if abs(g - GS[1][0]) < 1e-9 else None)
ax.axhline(math.sqrt(1.5), color="gray", ls="--", lw=1.0)
ax.text(0.02, math.sqrt(1.5) + 0.02, r"free $\sqrt{3/2}$", color="gray", fontsize=9,
        transform=ax.get_yaxis_transform())
ax.set_xlabel(r"$1/L^2$")
ax.set_ylabel(r"$\Delta_{F,\ell=2}/\Delta_{F,\ell=1}$")
ax.set_title(r"Glueball $F$ ratio $\ell=2/\ell=1$ vs $1/L^2$")
ax.set_xlim(-0.05, 1.1)
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("figs/glue_ratio_Fl2_Fl1_vs_invL2_claude.png", dpi=150)
print("wrote glue_ratio_Fl2_Fl1_vs_invL2_claude.png")
