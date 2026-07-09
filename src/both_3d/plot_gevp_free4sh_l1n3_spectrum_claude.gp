# plot_gevp_free4sh_l1n3_spectrum_claude.gp
# Free ensemble, face_sign ON, 4-shape basis {triangle, rect, fig8, three-triangle}, l=1 sector,
# rebased at trebase=1 -> 3 physical states (m=-1,0,+1 triplet), tmax=16, 10000 configs.
# cols: 1=t 2,3=ground(dup s2) then s=0..2 as (mean,err) -> cols 4..9.
#   gnuplot plot_gevp_free4sh_l1n3_spectrum_claude.gp -> gevp_free4sh_l1n3_spectrum_claude.png

set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set output "gevp_free4sh_l1n3_spectrum_claude.png"

set title "Free F_{12} l=1 spectrum -- 4-shape (+three-triangle), 3 states, 10k cfg"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.7]
set yrange [0.6:1.8]
set key top left

set arrow from 0.1,sqrt(2) to 2.7,sqrt(2) nohead dt 2 lc rgb "black"
set label "{/Symbol \326}2" at 2.5,sqrt(2)+0.04
set arrow from 0.1,1.33242 to 2.7,1.33242 nohead dt 3 lc rgb "dark-green"
set label "1.332 (L=1 det)" at 0.2,1.28 tc rgb "dark-green"

f = "gevp_msm_free4sh_l1_n3_claude.dat"
plot for [s=0:2] f u ($1+0.03*s-0.03):(column(4+2*s)):(column(5+2*s)) \
     w yerrorbars pt (s+6) ps 1.4 lc (s+1) t sprintf("state %d", s)
