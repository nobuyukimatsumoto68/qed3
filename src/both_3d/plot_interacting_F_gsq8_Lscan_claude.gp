# plot_interacting_F_gsq8_Lscan_claude.gp
# Headline: F (linear) l=1 effective mass vs t, gsq=8, L=1/2/4, one panel per Nf overlaid.
# Physical F l=1 = ground = col 2,3 (lightest of the l=1 triplet). Free anchor sqrt(2)=1.41421
# (L=1 lattice det 1.33242). Color+marker per L (colour-blind safe).
#   gnuplot plot_interacting_F_gsq8_Lscan_claude.gp -> gevp_interacting_F_gsq8_Lscan_claude.png

set terminal pngcairo size 1100,720 enhanced font "Helvetica,14"
set output "gevp_interacting_F_gsq8_Lscan_claude.png"

set title "F (linear) l=1 effective mass -- gsq=8, L scan (ground = lightest l=1)"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.3]
set yrange [0.0:2.0]
set key top right

set arrow from 0.1,sqrt(2) to 2.3,sqrt(2) nohead dt 2 lc rgb "black"
set label "{/Symbol \326}2 = 1.414 (free)" at 1.5,sqrt(2)+0.06
set arrow from 0.1,1.33242 to 2.3,1.33242 nohead dt 3 lc rgb "gray40"
set label "1.332 (L=1 det)" at 0.15,1.33242-0.07 tc rgb "gray40"

# Nf=2 red circle, Nf=4 blue square, Nf=6 dark-green triangle ; L by dx offset + point size
plot \
  "gevp_msm_Nf2_g8.000000_L1_claude.dat" u ($1-0.02):2:3 w yerrorbars pt 7 ps 1.2 lc rgb "red"        t "Nf2 L1", \
  "gevp_msm_Nf2_g8.000000_L2_claude.dat" u ($1-0.01):2:3 w yerrorbars pt 6 ps 1.2 lc rgb "red"        t "Nf2 L2", \
  "gevp_msm_Nf2_g8.000000_L4_claude.dat" u ($1+0.00):2:3 w yerrorbars pt 4 ps 1.2 lc rgb "red"        t "Nf2 L4", \
  "gevp_msm_Nf4_g8.000000_L1_claude.dat" u ($1+0.01):2:3 w yerrorbars pt 5 ps 1.2 lc rgb "blue"       t "Nf4 L1", \
  "gevp_msm_Nf4_g8.000000_L2_claude.dat" u ($1+0.02):2:3 w yerrorbars pt 5 ps 0.9 lc rgb "web-blue"   t "Nf4 L2", \
  "gevp_msm_Nf4_g8.000000_L4_claude.dat" u ($1+0.03):2:3 w yerrorbars pt 4 ps 0.9 lc rgb "web-blue"   t "Nf4 L4", \
  "gevp_msm_Nf6_g8.000000_L1_claude.dat" u ($1+0.04):2:3 w yerrorbars pt 9 ps 1.2 lc rgb "dark-green" t "Nf6 L1", \
  "gevp_msm_Nf6_g8.000000_L2_claude.dat" u ($1+0.05):2:3 w yerrorbars pt 11 ps 1.2 lc rgb "dark-green" t "Nf6 L2", \
  "gevp_msm_Nf6_g8.000000_L4_claude.dat" u ($1+0.06):2:3 w yerrorbars pt 10 ps 1.2 lc rgb "dark-green" t "Nf6 L4"
