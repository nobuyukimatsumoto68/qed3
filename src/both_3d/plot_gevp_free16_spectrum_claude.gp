# plot_gevp_free16_spectrum_claude.gp
# Free ensemble (non-compact Maxwell), face_sign ON, untwisted 3-shape basis (tri+rect+fig8),
# l=1..3, tmax=16 -> GEVP spectrum (8 eigenstate effmasses vs t) out to t=2.6.
# Data cols: 1=t 2,3=ground(dup s7) then s=0..7 as (mean,err) -> cols 4..19.
# References: l=1 -> sqrt(2)=1.414 (paper Delta_0); l=2 -> sqrt(6)=2.449.
#   gnuplot plot_gevp_free16_spectrum_claude.gp -> gevp_free16_spectrum_claude.png

set terminal pngcairo size 1000,680 enhanced font "Helvetica,14"
set output "gevp_free16_spectrum_claude.png"

set title "Free F_{12} GEVP spectrum -- ON, untwisted, l=1..3, tmax=16"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.7]
set yrange [0.8:3.2]
set key top right

set arrow from 0.1,sqrt(2) to 2.7,sqrt(2) nohead dt 2 lc rgb "black"
set label "{/Symbol \326}2" at 2.45,sqrt(2)+0.06
set arrow from 0.1,sqrt(6) to 2.7,sqrt(6) nohead dt 2 lc rgb "black"
set label "{/Symbol \326}6" at 2.45,sqrt(6)+0.06

f = "gevp_msm_free16_ON_full_claude.dat"
plot for [s=0:7] f u ($1+0.01*s-0.03):(column(4+2*s)):(column(5+2*s)) \
     w yerrorbars pt (s+1) ps 1.2 lc (s+1) t sprintf("state %d", s)
