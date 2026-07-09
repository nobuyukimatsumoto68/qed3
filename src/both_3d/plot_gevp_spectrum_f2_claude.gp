# plot_gevp_spectrum_f2_claude.gp
# F_12^2 (0++) glueball: full GEVP spectrum (all 8 eigenstate effmasses vs t),
# as in the original glue_ylms3 notebook (cell 76/111). One png per Nf.
# Data cols: 1=t 2,3=ground(dup of s7) then s=0..7 as (mean,err) pairs -> cols 4..19.
# NOTE: state 7 (the "ground", largest eigenvalue) is the residual constant/vacuum
#   mode (m -> 0); the physical 0++ sits in the excited states (lower eigenvalues).
# Production point: L=1, gsq=8.
#   gnuplot plot_gevp_spectrum_f2_claude.gp

set terminal pngcairo size 900,650 enhanced font "Helvetica,14"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:0.9]
set yrange [-0.2:3.5]
set key top right

do for [nf in "2 4 6"] {
  set output sprintf("gevp_f2_spectrum_Nf%s_L1gsq8_claude.png", nf)
  set title sprintf("F_{12}^2 (0++) GEVP spectrum (Nf=%s, L=1, gsq=8)", nf)
  f = sprintf("gevp_f2_Nf%s_gsq8.000000_L1_claude.dat", nf)
  plot for [s=0:7] f u ($1+0.006*s-0.02):(column(4+2*s)):(column(5+2*s)) \
       w yerrorbars pt (s+1) ps 1.2 lc (s+1) t sprintf("state %d", s)
}
