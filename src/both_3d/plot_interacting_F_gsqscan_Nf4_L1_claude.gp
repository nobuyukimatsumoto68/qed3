# auto-generated: F l=1 effmass vs t, gsq scan, Nf=4, L1.
set terminal pngcairo size 1050,700 enhanced font "Helvetica,14"
set output "gevp_interacting_F_gsqscan_Nf4_L1_claude.png"
set title "F l=1 (ground = lightest triplet) -- gsq scan, Nf=4, L=1"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.1]
set yrange [0.0:2.2]
set key top right
set arrow from 0.1,1.41421356 to 2.1,1.41421356 nohead dt 2 lc rgb "gray30"
set label "{/Symbol \326}2 (free)" at 1.3,1.41421356+0.06 tc rgb "gray30"
plot \
  "gevp_msm_Nf4_g1.000000_L1_claude.dat" u ($1-0.03):2:3 w yerrorbars pt 7 ps 1.1 lc rgb "dark-violet" t "gsq=1", \
  "gevp_msm_Nf4_g2.000000_L1_claude.dat" u ($1-0.02):2:3 w yerrorbars pt 5 ps 1.1 lc rgb "red" t "gsq=2", \
  "gevp_msm_Nf4_g2.200000_L1_claude.dat" u ($1-0.01):2:3 w yerrorbars pt 9 ps 1.1 lc rgb "orange" t "gsq=2.2", \
  "gevp_msm_Nf4_g2.500000_L1_claude.dat" u ($1+0.00):2:3 w yerrorbars pt 11 ps 1.1 lc rgb "gold" t "gsq=2.5", \
  "gevp_msm_Nf4_g4.000000_L1_claude.dat" u ($1+0.01):2:3 w yerrorbars pt 13 ps 1.1 lc rgb "web-green" t "gsq=4", \
  "gevp_msm_Nf4_g8.000000_L1_claude.dat" u ($1+0.02):2:3 w yerrorbars pt 4 ps 1.1 lc rgb "blue" t "gsq=8", \
  "gevp_msm_Nf4_g12.000000_L1_claude.dat" u ($1+0.03):2:3 w yerrorbars pt 6 ps 1.1 lc rgb "black" t "gsq=12"
