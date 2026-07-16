# plot_interacting_F_gsq8_Lscan_claude.gp
# F (linear) l=1 TRIPLET-AVERAGED effective mass vs t, gsq=8, L scan. Uses the unit-weight triplet
# average (gevp_msmAVG_*, cols 2,3), binsize=50. L4 DROPPED (pure noise; bs=50 undefined at Nc=68).
# Free-continuum anchor sqrt(2)=1.41421 (L=1 lattice det 1.33242). Colour per Nf, marker per L.
#   gnuplot plot_interacting_F_gsq8_Lscan_claude.gp -> gevp_interacting_F_gsq8_Lscan_claude.png

set terminal pngcairo size 1100,720 enhanced font "Helvetica,14"
set output "gevp_interacting_F_gsq8_Lscan_claude.png"

set title "F l=1 triplet-averaged -- gsq=8, L scan (L1, L2), bs=50"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.3]
set yrange [0.0:2.0]
set key top right

set arrow from 0.1,sqrt(2) to 2.3,sqrt(2) nohead dt 2 lc rgb "black"
set label "{/Symbol \326}2 = 1.414 (free, continuum)" at 1.2,sqrt(2)+0.06
set arrow from 0.1,1.33242 to 2.3,1.33242 nohead dt 3 lc rgb "gray40"
set label "1.332 (L=1 det)" at 0.15,1.33242-0.07 tc rgb "gray40"

# Nf2 red circle, Nf4 blue square, Nf6 dark-green triangle ; L by fill (filled=L1, open=L2) + dx
plot \
  "gevp_msmAVG_Nf2_g8.000000_L1_claude.dat" u ($1-0.02):2:3 w yerrorbars pt 7 ps 1.3 lc rgb "red"        t "Nf2 L1", \
  "gevp_msmAVG_Nf2_g8.000000_L2_claude.dat" u ($1-0.01):2:3 w yerrorbars pt 6 ps 1.3 lc rgb "red"        t "Nf2 L2", \
  "gevp_msmAVG_Nf4_g8.000000_L1_claude.dat" u ($1+0.00):2:3 w yerrorbars pt 5 ps 1.3 lc rgb "blue"       t "Nf4 L1", \
  "gevp_msmAVG_Nf4_g8.000000_L2_claude.dat" u ($1+0.01):2:3 w yerrorbars pt 4 ps 1.3 lc rgb "blue"       t "Nf4 L2", \
  "gevp_msmAVG_Nf6_g8.000000_L1_claude.dat" u ($1+0.02):2:3 w yerrorbars pt 9 ps 1.5 lc rgb "dark-green" t "Nf6 L1", \
  "gevp_msmAVG_Nf6_g8.000000_L2_claude.dat" u ($1+0.03):2:3 w yerrorbars pt 11 ps 1.5 lc rgb "dark-green" t "Nf6 L2"
