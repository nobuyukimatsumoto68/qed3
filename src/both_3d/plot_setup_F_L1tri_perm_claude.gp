# plot_setup_F_L1tri_perm_claude.gp
# F l=1, L=1, triangle only (NO GEVP), 3 m-states shown INDIVIDUALLY (no m-averaging), Nf2 g8.
# Diagonal per-(l=1,m) effmass; cols: t, m=-1(2,3), m=0(4,5), m=+1(6,7).
set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set output "gevp_setup_F_L1tri_perm_Nf2g8_claude.png"
set title "F l=1, L=1 -- triangle only (NO GEVP), 3 m-states (no averaging), Nf2 g8"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:1.3]
set yrange [1.05:1.55]
set key top right
set arrow 1 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2 (free, continuum)" at 0.72,sqrt(2)+0.01 tc rgb "gray30"
f = "reb_F_L1_tri_perm_Nf2g8_claude.dat"
plot f u ($1-0.01):2:3 w yerrorbars pt 7 ps 1.3 lc rgb "red"        t "m=-1", \
     f u ($1+0.00):4:5 w yerrorbars pt 5 ps 1.3 lc rgb "blue"       t "m=0", \
     f u ($1+0.01):6:7 w yerrorbars pt 9 ps 1.4 lc rgb "dark-green" t "m=+1"
