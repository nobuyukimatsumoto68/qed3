#!/usr/bin/env python3
# Per-m const-fit plots: data points + per-m fit band (M_s +/- sigma_s) with chi2/dof annotation +
# (for F) the inverse-variance averaged value as a shaded band. Fits via fit_perm_claude.
import subprocess
import numpy as np
from fit_perm_claude import load, fit_state

CH = [
  dict(name="F_L1", title="F l=1, L=1 -- per-m DIAG-chi2 fit + var-avg, Nf2 g8",
       eff="permF_L1_Nf2g8_claude.dat", jk="jkperm_F_L1_Nf2g8_claude.dat",
       states=[(0,"m=-1",'"red"',7),(1,"m=0",'"blue"',5),(2,"m=+1",'"dark-green"',9)],
       tlo=0.2, thi=0.8, mavg=True, anchor=("sqrt(2)",'{/Symbol \\326}2'),
       xr="0.1:1.3", yr="1.15:1.5"),
  dict(name="F_L2", title="F l=1, L=2 -- per-m DIAG-chi2 fit + var-avg (m=0 low), Nf2 g8",
       eff="permF_L2_Nf2g8_claude.dat", jk="jkperm_F_L2_Nf2g8_claude.dat",
       states=[(0,"m=-1",'"red"',7),(1,"m=0",'"blue"',5),(2,"m=+1",'"dark-green"',9)],
       tlo=0.2, thi=0.8, mavg=True, anchor=("sqrt(2)",'{/Symbol \\326}2'),
       xr="0.1:1.3", yr="0.8:1.7"),
  dict(name="F2_L1", title="F^2 0++, L=1 -- DIAG-chi2 fit, Nf2 g8",
       eff="permF2_L1_eff_Nf2g8_claude.dat", jk="jkperm_F2_L1_Nf2g8_claude.dat",
       states=[(0,"0++",'"blue"',7)],
       tlo=0.2, thi=0.6, mavg=False, anchor=("2*sqrt(2)",'2{/Symbol \\326}2'),
       xr="0.1:1.3", yr="-0.2:3.5"),
]

def band(oid, tlo, thi, lo, hi, color, alpha):
    return f'set object {oid} rect from {tlo},{lo} to {thi},{hi} fc rgb {color} fs transparent solid {alpha} noborder front'

for ch in CH:
    jk, at = load(ch["jk"])
    ntpts = jk.shape[1]
    ti_list = [ti for ti in range(ntpts) if ch["tlo"]-1e-9 <= (ti+1)*at <= ch["thi"]+1e-9]
    L = []
    L.append('set terminal pngcairo size 1050,680 enhanced font "Helvetica,14"')
    L.append(f'set output "gevp_fit_{ch["name"]}_Nf2g8_claude.png"')
    L.append(f'set title "{ch["title"]}"')
    L.append('set xlabel "t"')
    L.append('set ylabel "{/Symbol D}_{eff}"')
    L.append('set grid')
    L.append(f'set xrange [{ch["xr"]}]')
    L.append(f'set yrange [{ch["yr"]}]')
    L.append('set key top right')
    # free-limit anchor
    aval, alab = ch["anchor"]
    L.append(f'set arrow 100 from {ch["xr"].split(":")[0]},{aval} to {ch["xr"].split(":")[1]},{aval} nohead dt 2 lc rgb "gray40"')
    L.append(f'set label 100 "{alab}" at graph 0.82,first {aval} offset 0,0.6 tc rgb "gray40"')
    # per-state fits
    oid = 1
    Mbs, sigs = [], []
    for (s, lab, col, pt) in ch["states"]:
        Mb, M, sig, chi2, dof = fit_state(jk, ti_list, s)
        Mbs.append(Mb); sigs.append(sig)
        # fit band + line over [tlo,thi]
        L.append(band(oid, ch["tlo"], ch["thi"], M-sig, M+sig, col, 0.15)); oid += 1
        L.append(f'set arrow {200+s} from {ch["tlo"]},{M} to {ch["thi"]},{M} nohead lc rgb {col} lw 2')
        L.append(f'set label {300+s} "{lab}: {M:.3f}({sig:.3f})  {{/Symbol c}}^2/dof={chi2/dof:.2f}" '
                 f'at graph 0.03,graph {0.93-0.06*s} tc rgb {col}')
    # variance-averaged band (F only)
    if ch["mavg"]:
        Mbs = np.array(Mbs); nbins = jk.shape[0]
        u = np.array([1.0/sg**2 for sg in sigs]); u /= u.sum()
        Mc = (u[:,None]*Mbs).sum(0); M = Mc.mean()
        err = np.sqrt((nbins-1.0)/nbins*np.sum((Mc-M)**2))
        L.append(band(oid, ch["tlo"], ch["thi"], M-err, M+err, '"black"', 0.22)); oid += 1
        L.append(f'set arrow 250 from {ch["tlo"]},{M} to {ch["thi"]},{M} nohead lc rgb "black" lw 3 dt 1')
        L.append(f'set label 260 "var-avg = {M:.3f}({err:.3f})" at graph 0.03,graph 0.72 tc rgb "black"')
    # data
    plots = []
    for (s, lab, col, pt) in ch["states"]:
        dx = -0.012 + 0.012*s
        plots.append(f'  "{ch["eff"]}" u ($1{dx:+.3f}):{4+2*s}:{5+2*s} w yerrorbars pt {pt} ps 1.2 lc rgb {col} t "{lab}"')
    L.append("plot \\\n" + ", \\\n".join(plots))
    gp = f"plot_fit_{ch['name']}_claude.gp"
    open(gp,"w").write("\n".join(L)+"\n")
    subprocess.run(["gnuplot", gp], check=True)
    print("wrote", f"gevp_fit_{ch['name']}_Nf2g8_claude.png")
