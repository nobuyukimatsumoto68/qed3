# plot_interacting_L2_gsq8_claude.gp
# Dedicated L=2 effective masses at gsq=8, Nf=2/4/6 overlaid. Two PNGs (F and F^2), binsize=50.
# F l=1 = TRIPLET-AVERAGED (gevp_msmAVG_*, cols 2,3); F^2 0++ = state 0 (gevp_f2_*, cols 4,5).
# Free-continuum anchors: F -> sqrt(2)=1.41421 ; F^2 -> 2sqrt(2)=2.82843. Nconf: Nf2=1601, Nf4=518, Nf6=289.
# Colour+marker per Nf (colour-blind safe): Nf2 red circle, Nf4 blue square, Nf6 dark-green triangle.
#   gnuplot plot_interacting_L2_gsq8_claude.gp

set terminal pngcairo size 1050,700 enhanced font "Helvetica,14"
set grid
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set xrange [0.1:1.7]
set key top right

# ---- F (linear) l=1, triplet-averaged ----
set output "gevp_interacting_F_L2_gsq8_claude.png"
set title "F l=1 triplet-averaged -- L=2, gsq=8, bs=50"
set yrange [0.0:2.0]
set arrow 1 from 0.1,sqrt(2) to 1.7,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2 = 1.414 (free, continuum)" at 0.7,sqrt(2)+0.06 tc rgb "gray30"
plot \
  "gevp_msmAVG_Nf2_g8.000000_L2_claude.dat" u ($1-0.015):2:3 w yerrorbars pt 7 ps 1.3 lc rgb "red"        t "Nf2 (1601)", \
  "gevp_msmAVG_Nf4_g8.000000_L2_claude.dat" u ($1+0.000):2:3 w yerrorbars pt 5 ps 1.3 lc rgb "blue"       t "Nf4 (518)", \
  "gevp_msmAVG_Nf6_g8.000000_L2_claude.dat" u ($1+0.015):2:3 w yerrorbars pt 9 ps 1.5 lc rgb "dark-green" t "Nf6 (289)"

# ---- F^2 (0++) ----
set output "gevp_interacting_F2_L2_gsq8_claude.png"
set title "F^2 0++ effective mass -- L=2, gsq=8, bs=50 (physical, state 0)"
set yrange [-0.4:4.0]
set arrow 1 from 0.1,2*sqrt(2) to 1.7,2*sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "2{/Symbol \326}2 = 2.828 (free, continuum)" at 0.6,2*sqrt(2)+0.1 tc rgb "gray30"
plot \
  "gevp_f2_Nf2_g8.000000_L2_claude.dat" u ($1-0.015):4:5 w yerrorbars pt 7 ps 1.3 lc rgb "red"        t "Nf2 (1601)", \
  "gevp_f2_Nf4_g8.000000_L2_claude.dat" u ($1+0.000):4:5 w yerrorbars pt 5 ps 1.3 lc rgb "blue"       t "Nf4 (518)", \
  "gevp_f2_Nf6_g8.000000_L2_claude.dat" u ($1+0.015):4:5 w yerrorbars pt 9 ps 1.5 lc rgb "dark-green" t "Nf6 (289)"
