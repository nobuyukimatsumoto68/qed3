set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set output "gevp_f2_free4sh_l0_1st_spectrum_claude.png"
set title "Free F^2 (0++) l=0 -- no-vacsub, omit-const, LIGHTEST ONLY (1 state), 10k cfg"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:2.7]
set yrange [1.4:4.6]
set key top left
set arrow from 0.1,2*sqrt(2) to 2.7,2*sqrt(2) nohead dt 2 lc rgb "black"
set label "2{/Symbol \326}2 = 2.828 (two-photon 0++)" at 1.2,2*sqrt(2)-0.18
plot "gevp_f2_free4sh_l0_1st_claude.dat" u 1:2:3 w yerrorbars pt 7 ps 1.6 lc rgb "blue" t "physical 0++ (lightest only)"
