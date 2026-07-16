# plot_setup_p2_claude.gp
# Final phase-2 analysis setup, Nf2 g8: (1) F L=1 = triangle, NO GEVP, m-averaged; (2) F L=2 = GEVP
# rebase to 3 states; (3) F^2 L=1 = GEVP rebase to 2 (const/vacuum + ground 0++).
set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set key top right

# ---- (1) F L=1 triangle, m-averaged (no GEVP) ----
set output "gevp_setup_F_L1tri_mavg_Nf2g8_claude.png"
set title "F l=1, L=1 -- triangle only (NO GEVP), m-averaged, Nf2 g8"
set xrange [0.1:1.3]
set yrange [1.1:1.55]
set arrow 1 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2 (free, continuum)" at 0.75,sqrt(2)+0.01 tc rgb "gray30"
plot "reb_F_L1_tri_Nf2g8_claude.dat" u 1:2:3 w yerrorbars pt 7 ps 1.4 lc rgb "red" t "triangle m-avg"

# ---- (2) F L=2 full, rebase to 3 states ----
unset arrow 1; unset label 1
set output "gevp_setup_F_L2_3st_Nf2g8_claude.png"
set title "F l=1, L=2 full 5-shape -- GEVP rebase to 3 states, Nf2 g8"
set xrange [0.1:1.3]
set yrange [0.5:2.0]
set arrow 1 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2" at 1.15,sqrt(2)+0.04 tc rgb "gray30"
f2 = "reb_F_L2_3st_Nf2g8_claude.dat"
plot f2 u ($1-0.012):4:5 w yerrorbars pt 5 ps 1.2 lc rgb "blue"       t "state 0", \
     f2 u ($1+0.000):6:7 w yerrorbars pt 9 ps 1.3 lc rgb "dark-green" t "state 1", \
     f2 u ($1+0.012):8:9 w yerrorbars pt 7 ps 1.3 lc rgb "red"        t "ground (spurious)"

# ---- (3) F^2 L=1, rebase to 2 (const + ground 0++) ----
unset arrow 1; unset label 1
set output "gevp_setup_F2_L1_2st_Nf2g8_claude.png"
set title "F^2 0++, L=1 full 5-shape -- GEVP rebase to 2 (const + 0++), Nf2 g8"
set xrange [0.1:1.3]
set yrange [-0.4:3.5]
set arrow 1 from 0.1,2*sqrt(2) to 1.3,2*sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "2{/Symbol \326}2 (two-photon)" at 0.7,2*sqrt(2)+0.12 tc rgb "gray30"
set arrow 2 from 0.1,0 to 1.3,0 nohead dt 3 lc rgb "dark-green"
f3 = "reb_F2_L1_2st_Nf2g8_claude.dat"
plot f3 u ($1+0.008):4:5 w yerrorbars pt 7 ps 1.4 lc rgb "blue"     t "physical 0++", \
     f3 u ($1-0.008):2:3 w yerrorbars pt 6 ps 1.3 lc rgb "dark-red" t "const/vacuum mode"
