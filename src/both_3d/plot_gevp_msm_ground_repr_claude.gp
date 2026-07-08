# plot_gevp_msm_ground_repr_claude.gp
# Representative linear F_12 glueball: GEVP ground-state effmass vs t, Nf=2, L=1, gsq=8,
# WITH the restored ell=0 channel (n_lm=16), stride=1, ~3400 configs.
# cols: 1=t 2=ground mean 3=ground err.
#   gnuplot plot_gevp_msm_ground_repr_claude.gp -> gevp_msm_ground_Nf2_L1gsq8_claude.png

set terminal pngcairo size 900,650 enhanced font "Helvetica,14"
set output "gevp_msm_ground_Nf2_L1gsq8_claude.png"

set title "Linear F_{12} ground effmass -- Nf=2, L=1, gsq=8 (ell=0..3, stride=1)"
set xlabel "t"
set ylabel "a_t m_{eff}"
set key top left
set grid
set xrange [0.1:0.9]
set yrange [0.9:1.4]

plot "gevp_msm_Nf2_gsq8.000000_L1_claude.dat" u 1:2:3 w yerrorbars pt 7 ps 1.6 lc rgb "red" t "Nf=2 (ell=0..3)"
