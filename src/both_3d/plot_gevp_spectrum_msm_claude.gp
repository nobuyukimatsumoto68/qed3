# plot_gevp_spectrum_msm_claude.gp
# Linear F_12 glueball: full GEVP spectrum (all 8 eigenstate effmasses vs t),
# as in the original glue_ylms3 notebook (cell 76/111). One png per Nf.
# Data cols: 1=t 2,3=ground(dup of s7) then s=0..7 as (mean,err) pairs -> cols 4..19.
# Production point: L=1, gsq=8.
#   gnuplot plot_gevp_spectrum_msm_claude.gp

set terminal pngcairo size 900,650 enhanced font "Helvetica,14"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:0.9]
set yrange [0.8:2.3]
set key top right

do for [nf in "2 4 6"] {
  set output sprintf("gevp_msm_spectrum_Nf%s_L1gsq8_claude.png", nf)
  set title sprintf("Linear F_{12} GEVP spectrum (Nf=%s, L=1, gsq=8)", nf)
  f = sprintf("gevp_msm_Nf%s_gsq8.000000_L1_claude.dat", nf)
  plot for [s=0:7] f u ($1+0.006*s-0.02):(column(4+2*s)):(column(5+2*s)) \
       w yerrorbars pt (s+1) ps 1.2 lc (s+1) t sprintf("state %d", s)
}
