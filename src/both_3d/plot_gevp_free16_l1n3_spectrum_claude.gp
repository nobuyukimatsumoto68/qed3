# plot_gevp_free16_l1n3_spectrum_claude.gp
# Free ensemble, face_sign ON, l=1 sector (9 ops), rebased at trebase=1, reduced to the 3 PHYSICAL
# states (the m=-1,0,+1 triplet). tmax=16, out to t=2.6. cols: 1=t 2,3=ground(dup s2) then
# s=0..2 as (mean,err) -> cols 4..9. All 3 sit on sqrt(2)=1.414 (paper Delta_0; L=1 lattice ~1.33).
#   gnuplot plot_gevp_free16_l1n3_spectrum_claude.gp -> gevp_free16_l1n3_spectrum_claude.png

set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set output "gevp_free16_l1n3_spectrum_claude.png"

set title "Free F_{12} GEVP spectrum -- ON, l=1, 9->3 states (m-triplet), tmax=16"
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

f = "gevp_msm_free16_ON_l1_n3_claude.dat"
plot for [s=0:2] f u ($1+0.03*s-0.03):(column(4+2*s)):(column(5+2*s)) \
     w yerrorbars pt (s+6) ps 1.4 lc (s+1) t sprintf("state %d", s)
