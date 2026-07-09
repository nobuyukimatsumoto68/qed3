# plot_gevp_free16_l1n6_spectrum_claude.gp
# Free ensemble, face_sign ON, l=1 sector only (9 ops = 3 shapes x 3 m), rebased at trebase=1
# and reduced to 6 states, tmax=16. GEVP spectrum (6 eigenstate effmasses vs t) out to t=2.6.
# cols: 1=t 2,3=ground(dup s5) then s=0..5 as (mean,err) -> cols 4..15.
# l=1 primary -> sqrt(2)=1.414 (paper Delta_0); the m=-1,0,+1 triplet is 3-fold degenerate there.
#   gnuplot plot_gevp_free16_l1n6_spectrum_claude.gp -> gevp_free16_l1n6_spectrum_claude.png

set terminal pngcairo size 1000,680 enhanced font "Helvetica,14"
set output "gevp_free16_l1n6_spectrum_claude.png"

set title "Free F_{12} GEVP spectrum -- ON, l=1 only, 9->6 states, tmax=16"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.7]
set yrange [0.6:3.4]
set key top right

set arrow from 0.1,sqrt(2) to 2.7,sqrt(2) nohead dt 2 lc rgb "black"
set label "{/Symbol \326}2" at 2.5,sqrt(2)+0.08

f = "gevp_msm_free16_ON_l1_n6_claude.dat"
# wide x-offset so the (near-)degenerate levels don't hide under each other
plot for [s=0:5] f u ($1+0.035*s-0.088):(column(4+2*s)):(column(5+2*s)) \
     w yerrorbars pt (s+1) ps 1.1 lc (s+1) t sprintf("state %d", s)
