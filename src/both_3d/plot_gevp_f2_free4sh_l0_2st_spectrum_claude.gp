# plot_gevp_f2_free4sh_l0_2st_spectrum_claude.gp
# Free ensemble, F^2 (0++) l=0, 4-shape, NO vacuum subtraction, rebased to 2 states = the constant/
# vacuum mode + the lightest physical 0++. tmax=16, 10k cfg.
# cols: 1=t 2,3=ground(=constant mode, Delta~0) then s0 (=physical 0++) at cols 4,5.
# Physical 0++ = two-photon = 2*sqrt(2)=2.828.
#   gnuplot plot_gevp_f2_free4sh_l0_2st_spectrum_claude.gp -> gevp_f2_free4sh_l0_2st_spectrum_claude.png

set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set output "gevp_f2_free4sh_l0_2st_spectrum_claude.png"

set title "Free F^2 (0++) l=0 -- no-vacsub, 2 states (const + lightest), 10k cfg"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.7]
set yrange [-0.4:3.6]
set key center right

set arrow from 0.1,2*sqrt(2) to 2.7,2*sqrt(2) nohead dt 2 lc rgb "black"
set label "2{/Symbol \326}2 = 2.828 (two-photon 0++)" at 0.9,2*sqrt(2)+0.14
set arrow from 0.1,0 to 2.7,0 nohead dt 3 lc rgb "dark-green"
set label "0 (vacuum/constant mode)" at 0.2,-0.24 tc rgb "dark-green"

f = "gevp_f2_free4sh_l0_2st_claude.dat"
plot f u ($1+0.015):4:5 w yerrorbars pt 7 ps 1.6 lc rgb "blue"     t "physical 0++ (lightest)", \
     f u ($1-0.015):2:3 w yerrorbars pt 5 ps 1.4 lc rgb "dark-red" t "constant/vacuum mode"
