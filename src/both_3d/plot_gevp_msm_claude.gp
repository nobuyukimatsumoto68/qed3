# plot_gevp_msm_claude.gp
# Linear F_12 (pseudoscalar) glueball: GEVP ground-state effective mass vs t.
# Data cols: 1=t  2=ground mean  3=ground err  (ground = largest eigenvalue = lightest).
# Production point: gsq=8, L=1, Nf = 2 / 4 / 6.
# Colorblind-safe: distinct marker AND color per series.
#   gnuplot plot_gevp_msm_claude.gp   ->  gevp_msm_L1gsq8_claude.png

set terminal pngcairo size 900,650 enhanced font "Helvetica,14"
set output "gevp_msm_L1gsq8_claude.png"

set title "Linear F_{12} glueball -- GEVP ground effmass (L=1, gsq=8)"
set xlabel "t"
set ylabel "a_t m_{eff}"
set key top left
set grid
set xrange [0.1:0.9]
set yrange [0.8:1.4]

plot \
  "gevp_msm_Nf2_gsq8.000000_L1_claude.dat" u ($1-0.008):2:3 w yerrorbars pt 7 ps 1.4 lc rgb "red"       t "Nf=2", \
  "gevp_msm_Nf4_gsq8.000000_L1_claude.dat" u ($1+0.000):2:3 w yerrorbars pt 5 ps 1.4 lc rgb "blue"      t "Nf=4", \
  "gevp_msm_Nf6_gsq8.000000_L1_claude.dat" u ($1+0.008):2:3 w yerrorbars pt 9 ps 1.6 lc rgb "dark-green" t "Nf=6"
