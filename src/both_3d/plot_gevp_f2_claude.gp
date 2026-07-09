# plot_gevp_f2_claude.gp
# F_{12}^2 action-density (0++ scalar) glueball: GEVP ground-state effective mass vs t.
# Data cols: 1=t  2=ground mean  3=ground err  (ground = largest eigenvalue = lightest).
# NOTE: for F^2 the GEVP GROUND is the residual CONSTANT / vacuum mode (lambda -> 1),
#   so a_t m_eff -> 0. The physical 0++ is an EXCITED eigenvalue (not shown here).
#   This plot documents the vacuum contamination of the ground channel.
# Production point: gsq=8, L=1, Nf = 2 / 4 / 6.
# Colorblind-safe: distinct marker AND color per series.
#   gnuplot plot_gevp_f2_claude.gp   ->  gevp_f2_L1gsq8_claude.png

set terminal pngcairo size 900,650 enhanced font "Helvetica,14"
set output "gevp_f2_L1gsq8_claude.png"

set title "F_{12}^2 (0++) glueball -- GEVP ground effmass (L=1, gsq=8) [vacuum mode]"
set xlabel "t"
set ylabel "a_t m_{eff}"
set key top right
set grid
set xrange [0.1:0.9]
set yrange [-0.02:0.08]

plot \
  "gevp_f2_Nf2_gsq8.000000_L1_claude.dat" u ($1-0.008):2:3 w yerrorbars pt 7 ps 1.4 lc rgb "red"        t "Nf=2", \
  "gevp_f2_Nf4_gsq8.000000_L1_claude.dat" u ($1+0.000):2:3 w yerrorbars pt 5 ps 1.4 lc rgb "blue"       t "Nf=4", \
  "gevp_f2_Nf6_gsq8.000000_L1_claude.dat" u ($1+0.008):2:3 w yerrorbars pt 9 ps 1.6 lc rgb "dark-green" t "Nf=6"
