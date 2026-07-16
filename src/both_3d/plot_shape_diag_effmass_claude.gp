# plot_shape_diag_effmass_claude.gp
# Scrutiny of the linear-F rank-1 claim: per-shape DIAGONAL effective mass for triangle/rect/fig8,
# split by m (l=1), Nf4 g2 L1. COLOUR = m (the 3 curves that actually differ); MARKER = shape.
# The 3 shapes at fixed m land EXACTLY on top of each other (rank-1: shapes redundant); the 3
# m-bands are the only real structure. dx offset per shape so the coincident points are visible.
# cols: t, then (a=tri,rect,fig8) x (ilm=m-1,m0,m+1): mean,err = 2+2*(a*3+ilm), +1.
#   gnuplot plot_shape_diag_effmass_claude.gp -> gevp_shape_diag_effmass_claude.png

set terminal pngcairo size 1050,700 enhanced font "Helvetica,14"
set output "gevp_shape_diag_effmass_claude.png"
set title "Linear F: per-shape diagonal effmass by m (Nf4 g2 L1) -- shapes coincide, m splits"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:1.5]
set yrange [1.2:1.8]
set key top left
set arrow from 0.1,sqrt(2) to 1.5,sqrt(2) nohead dt 2 lc rgb "gray30"
set label "{/Symbol \326}2 (free, continuum)" at 1.0,sqrt(2)+0.02 tc rgb "gray30"

f = "shape_diag_effmass_Nf4g2L1_claude.dat"
# colour by m: m=-1 red, m=0 blue, m=+1 dark-green ; marker by shape: tri=circle rect=square fig8=triangle
plot \
  f u ($1-0.02):2:3   w yerrorbars pt 7  ps 1.2 lc rgb "red"        t "tri  m=-1", \
  f u ($1-0.02):4:5   w yerrorbars pt 7  ps 1.2 lc rgb "blue"       t "tri  m=0", \
  f u ($1-0.02):6:7   w yerrorbars pt 7  ps 1.2 lc rgb "dark-green" t "tri  m=+1", \
  f u ($1+0.00):8:9   w yerrorbars pt 5  ps 1.2 lc rgb "red"        t "rect m=-1", \
  f u ($1+0.00):10:11 w yerrorbars pt 5  ps 1.2 lc rgb "blue"       t "rect m=0", \
  f u ($1+0.00):12:13 w yerrorbars pt 5  ps 1.2 lc rgb "dark-green" t "rect m=+1", \
  f u ($1+0.02):14:15 w yerrorbars pt 9  ps 1.4 lc rgb "red"        t "fig8 m=-1", \
  f u ($1+0.02):16:17 w yerrorbars pt 9  ps 1.4 lc rgb "blue"       t "fig8 m=0", \
  f u ($1+0.02):18:19 w yerrorbars pt 9  ps 1.4 lc rgb "dark-green" t "fig8 m=+1"
