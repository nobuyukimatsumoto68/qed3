# plot_gevp_msm_spectrum_repr_claude.gp
# Representative linear F_12 glueball: full GEVP spectrum (8 eigenstate effmasses vs t),
# Nf=2, L=1, gsq=8, WITH the restored ell=0 channel. Style follows glue_ylms3 cell 76/111.
# cols: 1=t 2,3=ground(dup s7) then s=0..7 as (mean,err) -> cols 4..19.
#   gnuplot plot_gevp_msm_spectrum_repr_claude.gp -> gevp_msm_spectrum_Nf2_L1gsq8_claude.png

set terminal pngcairo size 900,650 enhanced font "Helvetica,14"
set output "gevp_msm_spectrum_Nf2_L1gsq8_claude.png"

set title "Linear F_{12} GEVP spectrum -- Nf=2, L=1, gsq=8 (ell=0..3)"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:0.9]
set yrange [0.8:2.3]
set key top right

f = "gevp_msm_Nf2_gsq8.000000_L1_claude.dat"
plot for [s=0:7] f u ($1+0.006*s-0.02):(column(4+2*s)):(column(5+2*s)) \
     w yerrorbars pt (s+1) ps 1.2 lc (s+1) t sprintf("state %d", s)
