# plot_gevp_msm_t0scan_claude.gp
# Linear F_12 ground effmass vs t, Nf=2 L=1 gsq=8 (ell=0..3), GEVP rebase t0 = trebase scanned
# over {0,1,2,3}. tcut=8 -> t = 0.2..1.4. cols: 1=t 2=ground mean 3=ground err.
# Colorblind-safe: distinct marker AND color per t0.
#   gnuplot plot_gevp_msm_t0scan_claude.gp -> gevp_msm_t0scan_Nf2_L1gsq8_claude.png

set terminal pngcairo size 1000,680 enhanced font "Helvetica,14"
set output "gevp_msm_t0scan_Nf2_L1gsq8_claude.png"

set title "Linear F_{12} ground effmass vs t0 (rebase) -- Nf=2, L=1, gsq=8 (ell=0..3)"
set xlabel "t"
set ylabel "a_t m_{eff}"
set key top left
set grid
set xrange [0.1:1.5]
set yrange [0.7:1.8]

plot \
  "gevp_msm_Nf2_gsq8_L1_t00_claude.dat" u ($1-0.012):2:3 w yerrorbars pt 7  ps 1.4 lc rgb "red"          t "t0=0", \
  "gevp_msm_Nf2_gsq8_L1_t01_claude.dat" u ($1-0.004):2:3 w yerrorbars pt 5  ps 1.4 lc rgb "blue"         t "t0=1", \
  "gevp_msm_Nf2_gsq8_L1_t02_claude.dat" u ($1+0.004):2:3 w yerrorbars pt 9  ps 1.6 lc rgb "dark-green"   t "t0=2", \
  "gevp_msm_Nf2_gsq8_L1_t03_claude.dat" u ($1+0.012):2:3 w yerrorbars pt 13 ps 1.6 lc rgb "dark-magenta" t "t0=3"
