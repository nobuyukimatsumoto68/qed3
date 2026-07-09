# plot_gevp_free16_l1_norebase_spectrum_claude.gp
# Free ensemble, face_sign ON, l=1 sector (9 ops), RAW full GEVP (NO rebase), all 9 states, tmax=16.
# cols: 1=t 2,3=ground(dup s8) then s=0..8 as (mean,err) -> cols 4..21.
# Only 3 states are physical (the l=1 sqrt2 triplet); the other 6 are the rank-3 null space
# (NaN / inf / huge Delta) and fall outside this y-range on purpose.
#   gnuplot plot_gevp_free16_l1_norebase_spectrum_claude.gp -> gevp_free16_l1_norebase_spectrum_claude.png

set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set output "gevp_free16_l1_norebase_spectrum_claude.png"

set title "Free F_{12} GEVP spectrum -- ON, l=1, RAW (no rebase), 9 states, tmax=16"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.7]
set yrange [0.6:1.9]
set key top left

set arrow from 0.1,sqrt(2) to 2.7,sqrt(2) nohead dt 2 lc rgb "black"
set label "{/Symbol \326}2" at 2.5,sqrt(2)+0.04
set arrow from 0.1,1.33242 to 2.7,1.33242 nohead dt 3 lc rgb "dark-green"
set label "1.332 (L=1 det)" at 0.2,1.28 tc rgb "dark-green"

f = "gevp_msm_free16_ON_l1_norebase_claude.dat"
plot for [s=0:8] f u ($1+0.02*s-0.08):(column(4+2*s)):(column(5+2*s)) \
     w yerrorbars pt (s+1) ps 1.2 lc (s+1) t sprintf("state %d", s)
