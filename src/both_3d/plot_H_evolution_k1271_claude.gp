#!/usr/bin/gnuplot
# Plot H evolution + Sf force norms + eta/phi ratio for normal trajectory k=1271.
# Data files:
#   H_evolution_k1271_claude.dat    -- cols: t type outer half inner_step Sg Sf_cached Smom H_approx
#   Sf_force_norms_k1271_claude.dat -- cols: t kick_idx pf_idx label L2_norm Linf_norm
#   eta_phi_k1271_claude.dat        -- cols: t kick_idx outer phase ipf eta_phi_ratio

set terminal pngcairo size 1400,2100 font "Sans,12"
set output "H_evolution_k1271_claude.png"

dat  = "H_evolution_k1271_claude.dat"
fdat = "Sf_force_norms_k1271_claude.dat"
edat = "eta_phi_k1271_claude.dat"

set multiplot layout 7,1 margins 0.10,0.97,0.06,0.98 spacing 0,0.03

# -----------------------------------------------------------------------
# Panel 1: H evolution
# -----------------------------------------------------------------------
set lmargin at screen 0.10
set rmargin at screen 0.97

set xlabel ""
set format x ""
set ylabel "H_{approx}"
set title "H evolution during normal trajectory k=1271  (dH = +0.567, rejected)"
set key top right

set xrange [-10:1015]
set yrange [10671:10682]

set style line 10 lc rgb "#aaaaaa" lw 1 dt 2
set for [t = 0:1010:101] arrow from t, graph 0 to t, graph 1 nohead ls 10

set label 1 "outer 0"  at  50, 10681.2 center font "Sans,10" tc rgb "#444444"
set label 2 "outer 1"  at 252, 10681.2 center font "Sans,10" tc rgb "#444444"
set label 3 "outer 2"  at 454, 10681.2 center font "Sans,10" tc rgb "#444444"
set label 4 "outer 3"  at 656, 10681.2 center font "Sans,10" tc rgb "#444444"
set label 5 "outer 4"  at 858, 10681.2 center font "Sans,10" tc rgb "#444444"
set label 6 "half 0|1" at 151, 10680.8 center font "Sans,9"  tc rgb "#888888"
set label 7 "half 0|1" at 353, 10680.8 center font "Sans,9"  tc rgb "#888888"
set label 8 "half 0|1" at 555, 10680.8 center font "Sans,9"  tc rgb "#888888"
set label 9 "half 0|1" at 757, 10680.8 center font "Sans,9"  tc rgb "#888888"
set label 10 "half 0|1" at 959, 10680.8 center font "Sans,9" tc rgb "#888888"

set style line 1 lc rgb "#888888" lw 1 pt 0
set style line 2 lc rgb "#0055cc" lw 0 pt 7 ps 1.2

set style line 7 lc rgb "#cc0000" lw 0 pt 5 ps 1.6

plot dat using ($2==0 ? $1 : 1/0):9 with lines  ls 1 title "H_{approx} (S_f cached)", \
     dat using ($2==1 ? $1 : 1/0):9 with points ls 2 title "PrintH (full S_f from CG)", \
     dat using ($2==3 ? $1 : 1/0):9 with points ls 7 title "h_0 (true, before kick)"

# -----------------------------------------------------------------------
# Panel 2: L2 norm (Frobenius) of Sf force at each kick
# -----------------------------------------------------------------------
unset label; unset arrow
set format x ""
set xlabel ""
set ylabel "||dSf||_2"
set title ""
set key top right

set yrange [9.0:11.5]
set xrange [-10:1015]

set style line 3 lc rgb "#0055cc" lw 1.5 pt 7 ps 1.0
set style line 4 lc rgb "#e07000" lw 1.5 pt 6 ps 1.0

set for [t = 0:1010:101] arrow from t, graph 0 to t, graph 1 nohead ls 10

plot fdat using ($3==0 ? $1 : 1/0):5 with linespoints ls 3 title "pf 0 L2", \
     fdat using ($3==1 ? $1 : 1/0):5 with linespoints ls 4 title "pf 1 L2"

# -----------------------------------------------------------------------
# Panel 3: Linf norm of Sf force at each kick
# -----------------------------------------------------------------------
unset label; unset arrow
set format x "%g"
set xlabel "step index t (kick positions at multiples of 101)"
set ylabel "||dSf||_{inf}"
set title ""
set key top right

set yrange [0.7:1.3]
set xrange [-10:1015]

set style line 5 lc rgb "#0055cc" lw 1.5 pt 7 ps 1.0
set style line 6 lc rgb "#e07000" lw 1.5 pt 6 ps 1.0

set for [t = 0:1010:101] arrow from t, graph 0 to t, graph 1 nohead ls 10

plot fdat using ($3==0 ? $1 : 1/0):6 with linespoints ls 5 title "pf 0 Linf", \
     fdat using ($3==1 ? $1 : 1/0):6 with linespoints ls 6 title "pf 1 Linf"

# -----------------------------------------------------------------------
# Panel 4: ||eta||/||phi|| ratio per pf at each kick
# -----------------------------------------------------------------------
unset label; unset arrow
set format x "%g"
set xlabel "step index t (kick positions at multiples of 101)"
set ylabel "||eta|| / ||phi||"
set title ""
set key top left

set yrange [0.26:0.32]
set xrange [-10:1015]

set style line 8  lc rgb "#0055cc" lw 1.5 pt 7 ps 1.0
set style line 9  lc rgb "#e07000" lw 1.5 pt 6 ps 1.0

set for [t = 0:1010:101] arrow from t, graph 0 to t, graph 1 nohead ls 10

plot edat using ($5==0 ? $1 : 1/0):6 with linespoints ls 8 title "pf 0", \
     edat using ($5==1 ? $1 : 1/0):6 with linespoints ls 9 title "pf 1"

# -----------------------------------------------------------------------
# Panel 5: Sg (gauge action)
# -----------------------------------------------------------------------
unset label; unset arrow
set format x ""
set xlabel ""
set ylabel "S_g"
set title ""
set key top right

set yrange [1450:2950]
set xrange [-10:1015]

set for [t = 0:1010:101] arrow from t, graph 0 to t, graph 1 nohead ls 10

plot dat using ($2==0 ? $1 : 1/0):6 with lines  ls 1 title "inner loop", \
     dat using ($2==1 ? $1 : 1/0):6 with points ls 2 title "PrintH"

# -----------------------------------------------------------------------
# Panel 6: Sf (fermion action, PrintH exact values only)
# -----------------------------------------------------------------------
unset label; unset arrow
set format x ""
set xlabel ""
set ylabel "S_f"
set title ""
set key top right

set yrange [6055:6125]
set xrange [-10:1015]

set for [t = 0:1010:101] arrow from t, graph 0 to t, graph 1 nohead ls 10

plot dat using ($2==1 ? $1 : 1/0):7 with points ls 2 title "PrintH (CG)", \
     dat using ($2==3 ? $1 : 1/0):7 with points ls 7 title "h_0 (true)"

# -----------------------------------------------------------------------
# Panel 7: Smom (momentum squared = (1/2) pi^2)
# -----------------------------------------------------------------------
unset label; unset arrow
set format x "%g"
set xlabel "step index t"
set ylabel "S_{mom}"
set title ""
set key top right

set yrange [1600:3100]
set xrange [-10:1015]

set for [t = 0:1010:101] arrow from t, graph 0 to t, graph 1 nohead ls 10

plot dat using ($2==0 ? $1 : 1/0):8 with lines  ls 1 title "inner loop", \
     dat using ($2==1 ? $1 : 1/0):8 with points ls 2 title "PrintH"

unset multiplot
