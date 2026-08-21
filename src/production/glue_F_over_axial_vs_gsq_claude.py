# glue_F_over_axial_vs_gsq_claude.py
# Glueball-F / axial-current mass ratios vs g^2:
#   R1 = Delta_{F,l=1} / Delta_{A,l=1}      R2 = Delta_{F,l=2} / Delta_{A,l=1}
# F l=1 and l=2 = GEVP per-m + fit_perm variance-avg (F l=1 [0.2,1.0] states 0,1,2; F l=2 (H) [0.2,0.6]
#   states 0..4).  Both GEVP effmasses are PHYSICAL (t-axis = a_t n_t; M ~ 1.3-1.7, same units as axial).
# axial l=1 = CANONICAL const+exp m0 read from axial_tp_masses_summary_claude.md (NOT the old const-fit
#   [4.0,5.2]) -- same definition used by the vector/axial ratio.  kmin=20.
# Free single-particle refs on S2: F l=1 = sqrt2, F l=2 = sqrt3, axial l=1 = 2 (Delta_l = l+1) ->
#   R1_free = sqrt2/2 = 1/sqrt2 ~ 0.707 ;  R2_free = sqrt3/2 ~ 0.866.
import subprocess, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NFS = [2, 4, 6]
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
Lmk = {1: "o", 2: "s", 3: "D", 4: "^"}
SUMMARY = "axial_tp_masses_summary_claude.md"


# ---------- axial l=1 (const+exp) from the summary ----------
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


# ---------- glueball mass via fit_perm ----------
def glue_mass(pref, tag, states, tlo, thi):
    jk = "%s_%s_jk_claude.dat" % (pref, tag)
    try:
        open(jk).close()
    except OSError:
        return None, None
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(tlo), str(thi), states],
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
            if not np.isfinite(v):
                return None, None
            return v, e
    return None, None


# Central window (PRIORITIZED) and alternate window; the shift central->alt is quoted as a SYSTEMATIC.
# central: F l=1 [0.2,1.0], F l=2 [0.2,0.6];  alt: both [0.4,1.4] (NM 2026-08-18).
CEN = {"F1": ("gevp_F", "0,1,2", 0.2, 1.0), "F2": ("gevp_Fl2", "0,1,2,3,4", 0.2, 0.6)}
ALT = {"F1": ("gevp_F", "0,1,2", 0.4, 1.4), "F2": ("gevp_Fl2", "0,1,2,3,4", 0.4, 1.4)}


# ---------- collect one ratio: central value, stat err, window systematic ----------
def collect(key):
    rows = []   # (L, nf, g, r, stat, sys, tot)
    for L in [1, 2, 3]:                 # L=4 dropped for F (NM 2026-08-18)
        for nf in NFS:
            for g in GS[L]:
                if (L, g, nf) not in axsum:
                    continue
                tag = "Nf%d_g%.6f_L%d" % (nf, g, L)
                mc, ec = glue_mass(CEN[key][0], tag, CEN[key][1], CEN[key][2], CEN[key][3])
                if mc is None:
                    continue
                a, ea = axsum[(L, g, nf)]
                r = mc / a
                stat = r * math.sqrt((ec / mc) ** 2 + (ea / a) ** 2)   # central-window statistical
                ma, ea2 = glue_mass(ALT[key][0], tag, ALT[key][1], ALT[key][2], ALT[key][3])
                sys = abs(r - ma / a) if ma is not None else 0.0        # window-variation systematic
                tot = math.sqrt(stat ** 2 + sys ** 2)
                rows.append((L, nf, g, r, stat, sys, tot))
    return rows


PANELS = [
    ("F1", collect("F1"), 1.0 / math.sqrt(2.0), r"free $1/\sqrt{2}$",
     r"$\Delta_{F,\ell=1}/\Delta_{A,\ell=1}$", "F l=1 / axial l=1"),
    ("F2", collect("F2"), math.sqrt(3.0) / 2.0, r"free $\sqrt{3}/2$",
     r"$\Delta_{F,\ell=2}/\Delta_{A,\ell=1}$", "F l=2 / axial l=1"),
]

# two-tier error bar: outer cap = total (stat + sys quadrature), inner = stat only.
fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
for ax, (key, rows, free, freelab, ylab, title) in zip(axs, PANELS):
    print("# %s : %d points  (value  stat  sys)" % (title, len(rows)))
    for L, nf, g, r, stat, sys, tot in rows:
        print("  L%d Nf%d g%.1f : %.4f (stat %.4f)(sys %.4f)" % (L, nf, g, r, stat, sys))
    for L in [1, 2, 3, 4]:
        for nf in NFS:
            d = sorted([x for x in rows if x[0] == L and x[1] == nf], key=lambda x: x[2])
            if not d:
                continue
            gg = [x[2] for x in d]
            rr = [x[3] for x in d]
            st = [x[4] for x in d]
            to = [x[6] for x in d]
            # outer total (thin, no marker)
            ax.errorbar(gg, rr, yerr=to, marker="none", color=nfcol[nf], capsize=4, ls="none",
                        lw=0.8, alpha=0.5)
            # inner stat (marker + line)
            ax.errorbar(gg, rr, yerr=st, marker=Lmk[L], color=nfcol[nf], capsize=2, ls="-",
                        lw=0.8, mfc="none" if L == 2 else nfcol[nf], label="Nf%d L%d" % (nf, L))
    ax.axhline(free, color="gray", ls="--", lw=1.0)
    ax.text(0.02, free, " " + freelab, color="gray", fontsize=9, va="bottom",
            transform=ax.get_yaxis_transform())
    ax.set_xlabel(r"$g_0^2$")
    ax.set_ylabel(ylab)
    ax.set_title(title + r" vs $g_0^2$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
fig.suptitle(r"Glueball $F$ ($\ell{=}1,2$) / axial $\ell{=}1$ vs $g_0^2$  (central = orig window; outer bar = +window systematic)")
fig.tight_layout()
fig.savefig("figs/ratio_F_over_axial_l12_vs_gsq_claude.png", dpi=150)
print("# wrote ratio_F_over_axial_l12_vs_gsq_claude.png")

# ---------- md: value(stat)(sys) ----------
out = ["# Glueball F(l=1,2)/axial(l=1) mass ratios vs g^2 (axial = const+exp summary)", "",
       "Central window: F l=1 [0.2,1.0], F l=2 [0.2,0.6].  Systematic = shift to [0.4,1.4].",
       "Format value(stat)(sys).", "",
       "| L | gsq | Nf | F1/A1 | F2/A1 |", "|---|-----|----|-------|-------|"]
d1 = {(L, nf, g): (r, s, y) for L, nf, g, r, s, y, t in PANELS[0][1]}
d2 = {(L, nf, g): (r, s, y) for L, nf, g, r, s, y, t in PANELS[1][1]}
for L, nf, g in sorted(set(d1) | set(d2)):
    s1 = "%.4f(%.4f)(%.4f)" % d1[(L, nf, g)] if (L, nf, g) in d1 else "--"
    s2 = "%.4f(%.4f)(%.4f)" % d2[(L, nf, g)] if (L, nf, g) in d2 else "--"
    out.append("| %d | %.1f | %d | %s | %s |" % (L, g, nf, s1, s2))
open("ratio_F_over_axial_l12_vs_gsq_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote ratio_F_over_axial_l12_vs_gsq_claude.md")
