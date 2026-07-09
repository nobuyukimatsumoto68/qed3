# plot_gevp_f2_free4sh_l0_spectrum_claude.gp
# Free ensemble, F^2 (0++) l=0 channel, 4-shape basis, vacuum-subtracted, tmax=16, 10k cfg.
# 4 GEVP states. cols: 1=t 2,3=ground(dup s3) then s=0..3 as (mean,err) -> cols 4..11.
# Expectation: physical 0++ = two-photon = 2*sqrt(2)=2.828 (state 2); a residual constant/vacuum
# mode survives at Delta~0 (the "ground", state 3); states 0,1 are null/noise directions.
#   gnuplot plot_gevp_f2_free4sh_l0_spectrum_claude.gp -> gevp_f2_free4sh_l0_spectrum_claude.png

set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set output "gevp_f2_free4sh_l0_spectrum_claude.png"

set title "Free F^2 (0++) l=0 spectrum -- 4-shape, vacuum-sub, 10k cfg"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.7]
set yrange [-0.4:3.6]
set key top right

set arrow from 0.1,2*sqrt(2) to 2.7,2*sqrt(2) nohead dt 2 lc rgb "black"
set label "2{/Symbol \326}2 (two-photon)" at 1.6,2*sqrt(2)+0.12
set arrow from 0.1,0 to 2.7,0 nohead dt 3 lc rgb "dark-green"
set label "0 (vacuum mode)" at 0.2,-0.22 tc rgb "dark-green"

f = "gevp_f2_free4sh_l0_claude.dat"
plot for [s=0:3] f u ($1+0.03*s-0.045):(column(4+2*s)):(column(5+2*s)) \
     w yerrorbars pt (s+6) ps 1.4 lc (s+1) t sprintf("state %d", s)
