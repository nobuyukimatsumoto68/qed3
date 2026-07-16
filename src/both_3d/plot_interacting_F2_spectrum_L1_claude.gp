# plot_interacting_F2_spectrum_L1_claude.gp
# F^2 (0++) GEVP spectrum at L=1, one PNG per Nf, for the reliable couplings gsq2 & gsq8 (bs=50).
# The l=0 F^2 GEVP keeps 2 states: the physical 0++ (state 0, cols 4,5 -> two-photon 2sqrt2) and the
# constant/vacuum mode (ground, cols 2,3 -> Delta~0, which the GEVP separates on its own).
# Filled = physical 0++ ; open = constant/vacuum mode. Colour per coupling (red gsq2, blue gsq8).
#   gnuplot plot_interacting_F2_spectrum_L1_claude.gp   (emits Nf2/Nf4/Nf6 PNGs)

set terminal pngcairo size 1050,700 enhanced font "Helvetica,14"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.1]
set yrange [-0.4:4.0]
set key top right

do for [nf in "2 4 6"] {
  set output sprintf("gevp_interacting_F2_spectrum_Nf%s_L1_claude.png", nf)
  set title sprintf("F^2 0++ GEVP spectrum -- L=1, Nf=%s, bs=50 (gsq2 & gsq8)", nf)
  set arrow 1 from 0.1,2*sqrt(2) to 2.1,2*sqrt(2) nohead dt 2 lc rgb "gray30"
  set label 1 "2{/Symbol \326}2 = 2.828 (free, continuum)" at 1.05,2*sqrt(2)+0.12 tc rgb "gray30"
  set arrow 2 from 0.1,0 to 2.1,0 nohead dt 3 lc rgb "dark-green"
  set label 2 "0 (vacuum/constant mode)" at 0.2,-0.25 tc rgb "dark-green"
  f2 = sprintf("gevp_f2_Nf%s_g2.000000_L1_claude.dat", nf)
  f8 = sprintf("gevp_f2_Nf%s_g8.000000_L1_claude.dat", nf)
  plot \
    f2 u ($1-0.02):4:5 w yerrorbars pt 7 ps 1.3 lc rgb "red"  t "0++ (gsq2)", \
    f2 u ($1-0.02):2:3 w yerrorbars pt 6 ps 1.3 lc rgb "red"  t "const (gsq2)", \
    f8 u ($1+0.02):4:5 w yerrorbars pt 5 ps 1.3 lc rgb "blue" t "0++ (gsq8)", \
    f8 u ($1+0.02):2:3 w yerrorbars pt 4 ps 1.3 lc rgb "blue" t "const (gsq8)"
}
