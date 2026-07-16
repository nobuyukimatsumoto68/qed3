# plot_reb_p2_claude.gp
# Rebased F l=1 GEVP, Nf2 g8. (1) L=2 full, rebase to 3 states; (2) L=1 full, rebase to 2 states.
# cols: t, ground(2,3), then s0,s1,... at (4,5),(6,7),...  (last s = ground = lightest).
set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set key top right

# ---- (1) F L=2 full, rebase to 3 states ----
set output "gevp_reb_F_L2_3st_Nf2g8_claude.png"
set title "F l=1, L=2 full 5-shape -- rebase to 3 states, Nf2 g8"
set xrange [0.1:1.3]
set yrange [0.5:2.0]
set arrow 1 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2" at 1.15,sqrt(2)+0.04 tc rgb "gray30"
f2 = "reb_F_L2_3st_Nf2g8_claude.dat"
plot f2 u ($1-0.012):4:5 w yerrorbars pt 5  ps 1.2 lc rgb "blue"       t "state 0", \
     f2 u ($1+0.000):6:7 w yerrorbars pt 9  ps 1.3 lc rgb "dark-green" t "state 1", \
     f2 u ($1+0.012):8:9 w yerrorbars pt 7  ps 1.3 lc rgb "red"        t "ground (lightest = spurious)"

# ---- (2) F L=1 full, rebase to 2 states ----
unset arrow 1; unset label 1
set output "gevp_reb_F_L1_2st_Nf2g8_claude.png"
set title "F l=1, L=1 full 5-shape -- rebase to 2 states, Nf2 g8"
set xrange [0.1:1.3]
set yrange [1.0:1.6]
set arrow 1 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2 (free, continuum)" at 0.75,sqrt(2)+0.02 tc rgb "gray30"
f1 = "reb_F_L1_2st_Nf2g8_claude.dat"
plot f1 u ($1-0.008):4:5 w yerrorbars pt 5 ps 1.2 lc rgb "blue" t "state 0 (stable ~1.36)", \
     f1 u ($1+0.008):6:7 w yerrorbars pt 7 ps 1.3 lc rgb "red"  t "ground (drifts down)"
