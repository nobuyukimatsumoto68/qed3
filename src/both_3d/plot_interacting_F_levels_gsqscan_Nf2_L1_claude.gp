set terminal pngcairo size 1050,700 enhanced font "Helvetica,14"
set output "gevp_interacting_F_levels_gsqscan_Nf2_L1_claude.png"
set title "F l=1 -- three levels (m-triplet) per coupling, Nf=2, L=1, bs=50"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.1]
set yrange [0.0:2.2]
set key top right
set arrow from 0.1,sqrt(2) to 2.1,sqrt(2) nohead dt 2 lc rgb "gray30"
set label "{/Symbol \326}2 = 1.414 (free, continuum)" at 1.15,sqrt(2)+0.06 tc rgb "gray30"
plot \
  "gevp_msm_Nf2_g1.000000_L1_claude.dat" u ($1-0.03):4:5 w yerrorbars pt 7 ps 1.0 lc rgb "dark-violet" t "gsq=1", \
  "gevp_msm_Nf2_g1.000000_L1_claude.dat" u ($1-0.03):6:7 w yerrorbars pt 7 ps 1.0 lc rgb "dark-violet" notitle, \
  "gevp_msm_Nf2_g1.000000_L1_claude.dat" u ($1-0.03):8:9 w yerrorbars pt 7 ps 1.0 lc rgb "dark-violet" notitle, \
  "gevp_msm_Nf2_g2.000000_L1_claude.dat" u ($1-0.02):4:5 w yerrorbars pt 5 ps 1.0 lc rgb "red" t "gsq=2", \
  "gevp_msm_Nf2_g2.000000_L1_claude.dat" u ($1-0.02):6:7 w yerrorbars pt 5 ps 1.0 lc rgb "red" notitle, \
  "gevp_msm_Nf2_g2.000000_L1_claude.dat" u ($1-0.02):8:9 w yerrorbars pt 5 ps 1.0 lc rgb "red" notitle, \
  "gevp_msm_Nf2_g2.400000_L1_claude.dat" u ($1-0.01):4:5 w yerrorbars pt 9 ps 1.0 lc rgb "orange" t "gsq=2.4", \
  "gevp_msm_Nf2_g2.400000_L1_claude.dat" u ($1-0.01):6:7 w yerrorbars pt 9 ps 1.0 lc rgb "orange" notitle, \
  "gevp_msm_Nf2_g2.400000_L1_claude.dat" u ($1-0.01):8:9 w yerrorbars pt 9 ps 1.0 lc rgb "orange" notitle, \
  "gevp_msm_Nf2_g2.500000_L1_claude.dat" u ($1+0.00):4:5 w yerrorbars pt 11 ps 1.0 lc rgb "gold" t "gsq=2.5", \
  "gevp_msm_Nf2_g2.500000_L1_claude.dat" u ($1+0.00):6:7 w yerrorbars pt 11 ps 1.0 lc rgb "gold" notitle, \
  "gevp_msm_Nf2_g2.500000_L1_claude.dat" u ($1+0.00):8:9 w yerrorbars pt 11 ps 1.0 lc rgb "gold" notitle, \
  "gevp_msm_Nf2_g4.000000_L1_claude.dat" u ($1+0.01):4:5 w yerrorbars pt 13 ps 1.0 lc rgb "web-green" t "gsq=4", \
  "gevp_msm_Nf2_g4.000000_L1_claude.dat" u ($1+0.01):6:7 w yerrorbars pt 13 ps 1.0 lc rgb "web-green" notitle, \
  "gevp_msm_Nf2_g4.000000_L1_claude.dat" u ($1+0.01):8:9 w yerrorbars pt 13 ps 1.0 lc rgb "web-green" notitle, \
  "gevp_msm_Nf2_g8.000000_L1_claude.dat" u ($1+0.02):4:5 w yerrorbars pt 4 ps 1.0 lc rgb "blue" t "gsq=8", \
  "gevp_msm_Nf2_g8.000000_L1_claude.dat" u ($1+0.02):6:7 w yerrorbars pt 4 ps 1.0 lc rgb "blue" notitle, \
  "gevp_msm_Nf2_g8.000000_L1_claude.dat" u ($1+0.02):8:9 w yerrorbars pt 4 ps 1.0 lc rgb "blue" notitle, \
  "gevp_msm_Nf2_g12.000000_L1_claude.dat" u ($1+0.03):4:5 w yerrorbars pt 6 ps 1.0 lc rgb "black" t "gsq=12", \
  "gevp_msm_Nf2_g12.000000_L1_claude.dat" u ($1+0.03):6:7 w yerrorbars pt 6 ps 1.0 lc rgb "black" notitle, \
  "gevp_msm_Nf2_g12.000000_L1_claude.dat" u ($1+0.03):8:9 w yerrorbars pt 6 ps 1.0 lc rgb "black" notitle
