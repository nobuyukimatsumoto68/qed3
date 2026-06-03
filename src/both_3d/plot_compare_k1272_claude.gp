#!/usr/bin/gnuplot
# Comparison of k=1272 trajectory: nsteps=5 (dH=-9.25) vs nsteps=7 (dH=-0.25) vs nsteps=10 (dH=-0.26).
# X axis: physical trajectory time t in [0, tmax=1.9].
# Data files:
#   H_compare_k1272_nsteps5_claude.dat  -- cols: t_phys ph_step Sg Sf Smom H
#   H_compare_k1272_nsteps7_claude.dat  -- same
#   H_compare_k1272_nsteps10_claude.dat -- same

set terminal pngcairo size 1200,1200 font "Sans,13"
set output "compare_k1272_claude.png"

dat5  = "H_compare_k1272_nsteps5_claude.dat"
dat7  = "H_compare_k1272_nsteps7_claude.dat"
dat10 = "H_compare_k1272_nsteps10_claude.dat"

set style line 1 lc rgb "#0055cc" lw 2 pt 7 ps 1.1   # nsteps=5  (blue)
set style line 2 lc rgb "#009933" lw 2 pt 4 ps 1.1   # nsteps=7  (green)
set style line 3 lc rgb "#cc3300" lw 2 pt 6 ps 1.1   # nsteps=10 (red)

set multiplot layout 4,1 margins 0.11,0.97,0.07,0.97 spacing 0,0.03

# -----------------------------------------------------------------------
# Panel 1: H
# -----------------------------------------------------------------------
set xlabel ""
set format x ""
set ylabel "H"
set title "k=1272: nsteps=5 (dH=-9.25) vs nsteps=7 (dH=-0.25) vs nsteps=10 (dH=-0.26)"
set key top right

set xrange [-0.05:1.95]
set yrange [10948:10966]

set label 1 "nsteps=5,  dH=-9.25"  at 1.85, 10951.5 right font "Sans,12" tc rgb "#0055cc"
set label 2 "nsteps=7,  dH=-0.25"  at 1.85, 10957.0 right font "Sans,12" tc rgb "#009933"
set label 3 "nsteps=10, dH=-0.26"  at 1.85, 10963.5 right font "Sans,12" tc rgb "#cc3300"

# h0 shared horizontal reference
set arrow 1 from 0, 10959.37 to 1.9, 10959.37 nohead lw 1 dt 3 lc rgb "#888888"
set label 4 "h_0" at 0.02, 10959.6 left font "Sans,11" tc rgb "#888888"

plot dat5  using 1:6 with linespoints ls 1 title "nsteps=5", \
     dat7  using 1:6 with linespoints ls 2 title "nsteps=7", \
     dat10 using 1:6 with linespoints ls 3 title "nsteps=10"

# -----------------------------------------------------------------------
# Panel 2: Sg
# -----------------------------------------------------------------------
unset label; unset arrow
set format x ""
set xlabel ""
set ylabel "S_g"
set title ""
set key top right

set xrange [-0.05:1.95]
set yrange [1930:2070]

plot dat5  using 1:3 with linespoints ls 1 title "nsteps=5", \
     dat7  using 1:3 with linespoints ls 2 title "nsteps=7", \
     dat10 using 1:3 with linespoints ls 3 title "nsteps=10"

# -----------------------------------------------------------------------
# Panel 3: Sf
# -----------------------------------------------------------------------
unset label; unset arrow
set format x ""
set xlabel ""
set ylabel "S_f"
set title ""
set key top right

set xrange [-0.05:1.95]
set yrange [6308:6332]

plot dat5  using 1:4 with linespoints ls 1 title "nsteps=5", \
     dat7  using 1:4 with linespoints ls 2 title "nsteps=7", \
     dat10 using 1:4 with linespoints ls 3 title "nsteps=10"

# -----------------------------------------------------------------------
# Panel 4: Smom
# -----------------------------------------------------------------------
unset label; unset arrow
set format x "%g"
set xlabel "trajectory time t"
set ylabel "S_{mom}"
set title ""
set key top right

set xrange [-0.05:1.95]
set yrange [2565:2690]

plot dat5  using 1:5 with linespoints ls 1 title "nsteps=5", \
     dat7  using 1:5 with linespoints ls 2 title "nsteps=7", \
     dat10 using 1:5 with linespoints ls 3 title "nsteps=10"

unset multiplot
