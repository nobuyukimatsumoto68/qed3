#!/usr/bin/env python3
# Generate gnuplot scripts for the interacting shape-basis GEVP effective masses and render them.
#   F (linear) l=1: per-state dat cols  ground(2,3) s0(4,5) s1(6,7) s2(8,9)  [3 kept states = m-triplet]
#   F l=1 averaged: gevp_msmAVG_* dat, unit-weight triplet mean per jk sample -> cols (2,3)
#   F^2 (0++)     : physical state 0 = cols (4,5)  [l=0, non-degenerate]
# Colour+marker per gsq (colour-blind safe: distinct pt AND lc per coupling).
# Anchors labelled "free, continuum": sqrt(2) / 2sqrt(2) are the CONTINUUM free-limit values
# (the L=1 lattice free value is 1.332, distinct).
import subprocess

GSTYLE = {
    "1.000000":  ('"dark-violet"', 7),
    "2.000000":  ('"red"',         5),
    "2.200000":  ('"orange"',      9),
    "2.400000":  ('"orange"',      9),
    "2.500000":  ('"gold"',        11),
    "4.000000":  ('"web-green"',   13),
    "8.000000":  ('"blue"',        4),
    "12.000000": ('"black"',       6),
}
GLABEL = {k: k.rstrip("0").rstrip(".") for k in GSTYLE}
NF_GSQ = {
    2: ["1.000000","2.000000","2.400000","2.500000","4.000000","8.000000","12.000000"],
    4: ["1.000000","2.000000","2.200000","2.500000","4.000000","8.000000","12.000000"],
    6: ["1.000000","2.000000","2.400000","4.000000","8.000000","12.000000"],
}

def run(fname):
    subprocess.run(["gnuplot", fname], check=True)

def header(png, title, ylab, yrange, xrange="0.1:2.1"):
    return [
        'set terminal pngcairo size 1050,700 enhanced font "Helvetica,14"',
        f'set output "{png}"',
        f'set title "{title}"',
        'set xlabel "t"',
        f'set ylabel "{ylab}"',
        'set grid',
        f'set xrange [{xrange}]',
        f'set yrange [{yrange}]',
        'set key top right',
    ]

def anchor(val, label):
    return [
        f'set arrow from 0.1,{val} to 2.1,{val} nohead dt 2 lc rgb "gray30"',
        f'set label "{label}" at 1.15,{val}+0.06 tc rgb "gray30"',
    ]

# ---- F l=1: THREE LEVELS (s0,s1,s2 same colour+marker per coupling) --------
for nf in (2,4,6):
    png = f"gevp_interacting_F_levels_gsqscan_Nf{nf}_L1_claude.png"
    L = header(png, f"F l=1 -- three levels (m-triplet) per coupling, Nf={nf}, L=1, bs=50",
               "{/Symbol D}_{eff}", "0.0:2.2")
    L += anchor("sqrt(2)", "{/Symbol \\326}2 = 1.414 (free, continuum)")
    plots = []
    for i,g in enumerate(NF_GSQ[nf]):
        lc,pt = GSTYLE[g]
        dx = -0.03 + 0.01*i
        dat = f"gevp_msm_Nf{nf}_g{g}_L1_claude.dat"
        # three states, identical colour+marker; only s0 carries the legend title
        plots.append(f'  "{dat}" u ($1{dx:+.2f}):4:5 w yerrorbars pt {pt} ps 1.0 lc rgb {lc} t "gsq={GLABEL[g]}"')
        plots.append(f'  "{dat}" u ($1{dx:+.2f}):6:7 w yerrorbars pt {pt} ps 1.0 lc rgb {lc} notitle')
        plots.append(f'  "{dat}" u ($1{dx:+.2f}):8:9 w yerrorbars pt {pt} ps 1.0 lc rgb {lc} notitle')
    L.append("plot \\\n" + ", \\\n".join(plots))
    gp = f"plot_interacting_F_levels_gsqscan_Nf{nf}_L1_claude.gp"
    open(gp,"w").write("\n".join(L)+"\n")
    run(gp); print("wrote", png)

# ---- F l=1: AVERAGED triplet (unit-weight mean per jk sample) ---------------
for nf in (2,4,6):
    png = f"gevp_interacting_F_avg_gsqscan_Nf{nf}_L1_claude.png"
    L = header(png, f"F l=1 -- unit-weight triplet average, Nf={nf}, L=1, bs=50",
               "{/Symbol D}_{eff}", "0.0:2.2")
    L += anchor("sqrt(2)", "{/Symbol \\326}2 = 1.414 (free, continuum)")
    plots = []
    for i,g in enumerate(NF_GSQ[nf]):
        lc,pt = GSTYLE[g]
        dx = -0.03 + 0.01*i
        dat = f"gevp_msmAVG_Nf{nf}_g{g}_L1_claude.dat"
        plots.append(f'  "{dat}" u ($1{dx:+.2f}):2:3 w yerrorbars pt {pt} ps 1.1 lc rgb {lc} t "gsq={GLABEL[g]}"')
    L.append("plot \\\n" + ", \\\n".join(plots))
    gp = f"plot_interacting_F_avg_gsqscan_Nf{nf}_L1_claude.gp"
    open(gp,"w").write("\n".join(L)+"\n")
    run(gp); print("wrote", png)

# ---- F l=1: INVERSE-VARIANCE weighted triplet average ----------------------
for nf in (2,4,6):
    png = f"gevp_interacting_F_wavg_gsqscan_Nf{nf}_L1_claude.png"
    L = header(png, f"F l=1 -- inverse-variance weighted triplet average, Nf={nf}, L=1, bs=50",
               "{/Symbol D}_{eff}", "0.0:2.2")
    L += anchor("sqrt(2)", "{/Symbol \\326}2 = 1.414 (free, continuum)")
    plots = []
    for i,g in enumerate(NF_GSQ[nf]):
        lc,pt = GSTYLE[g]
        dx = -0.03 + 0.01*i
        dat = f"gevp_msmWAVG_Nf{nf}_g{g}_L1_claude.dat"
        plots.append(f'  "{dat}" u ($1{dx:+.2f}):2:3 w yerrorbars pt {pt} ps 1.1 lc rgb {lc} t "gsq={GLABEL[g]}"')
    L.append("plot \\\n" + ", \\\n".join(plots))
    gp = f"plot_interacting_F_wavg_gsqscan_Nf{nf}_L1_claude.gp"
    open(gp,"w").write("\n".join(L)+"\n")
    run(gp); print("wrote", png)

# ---- F^2 0++ (l=0, non-degenerate): single physical state ------------------
for nf in (2,4,6):
    png = f"gevp_interacting_F2_gsqscan_Nf{nf}_L1_claude.png"
    L = header(png, f"F^2 0++ (physical, state 0) -- gsq scan, Nf={nf}, L=1, bs=50",
               "{/Symbol D}_{eff}", "-0.2:4.0")
    L += anchor("2*sqrt(2)", "2{/Symbol \\326}2 = 2.828 (free, continuum)")
    plots = []
    for i,g in enumerate(NF_GSQ[nf]):
        lc,pt = GSTYLE[g]
        dx = -0.03 + 0.01*i
        dat = f"gevp_f2_Nf{nf}_g{g}_L1_claude.dat"
        plots.append(f'  "{dat}" u ($1{dx:+.2f}):4:5 w yerrorbars pt {pt} ps 1.1 lc rgb {lc} t "gsq={GLABEL[g]}"')
    L.append("plot \\\n" + ", \\\n".join(plots))
    gp = f"plot_interacting_F2_gsqscan_Nf{nf}_L1_claude.gp"
    open(gp,"w").write("\n".join(L)+"\n")
    run(gp); print("wrote", png)

print("ALL PLOTS DONE")
