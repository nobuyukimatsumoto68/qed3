# plot_binscan_claude.gp
# Jackknife error on the averaged F l=1 effmass vs binsize, for the two highest-stat ensembles
# (Nf4 g2 Nc=5899, Nf6 g2 Nc=5538). The error saturates once bins decorrelate (autocorrelation);
# the saturated value is the honest error. Colour+marker per ensemble; open/filled per t.
#   gnuplot plot_binscan_claude.gp -> gevp_binscan_claude.png

set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set output "gevp_binscan_claude.png"
set title "Jackknife error vs binsize -- averaged F l=1 (autocorrelation saturation)"
set xlabel "binsize (configs per bin)"
set ylabel "jackknife error on {/Symbol D}_{eff}"
set grid
set logscale x
set xrange [0.8:320]
set yrange [0.004:0.023]
set key bottom right

set arrow from 50,0.004 to 50,0.023 nohead dt 2 lc rgb "gray50"
set label "chosen bs=50" at 55,0.0055 tc rgb "gray30"

plot \
  "binscan_Nf4g2_claude.dat" u 1:2 w linespoints pt 7  ps 1.3 lc rgb "red"        t "Nf4 g2, t=0.2", \
  "binscan_Nf4g2_claude.dat" u 1:3 w linespoints pt 6  ps 1.3 lc rgb "red"        t "Nf4 g2, t=0.4", \
  "binscan_Nf6g2_claude.dat" u 1:2 w linespoints pt 5  ps 1.3 lc rgb "dark-green" t "Nf6 g2, t=0.2", \
  "binscan_Nf6g2_claude.dat" u 1:3 w linespoints pt 4  ps 1.3 lc rgb "dark-green" t "Nf6 g2, t=0.4"
