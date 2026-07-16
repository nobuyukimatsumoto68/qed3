# plot_perm_F_claude.gp
# Explicit per-m GEVP for linear F l=1, ground of each m-block, Nf2 g8. m sorted ascending: -1,0,+1
# -> cols s0(4,5)=m=-1, s1(6,7)=m=0, s2(8,9)=m=+1.
#  (1) L=1: 3 m-grounds ~degenerate (statistical spread).
#  (2) L=2: m=-1,+1 agree ~1.31; m=0 alone is the spurious low mode -> breaks degeneracy.
set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set key top right

set output "gevp_perm_F_L1_Nf2g8_claude.png"
set title "F l=1, L=1 -- explicit per-m GEVP (ground of each m-block), Nf2 g8"
set xrange [0.1:1.3]
set yrange [1.05:1.55]
set arrow 1 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2" at 1.15,sqrt(2)+0.02 tc rgb "gray30"
f1 = "permF_L1_Nf2g8_claude.dat"
plot f1 u ($1-0.01):4:5 w yerrorbars pt 7 ps 1.3 lc rgb "red"        t "m=-1", \
     f1 u ($1+0.00):6:7 w yerrorbars pt 5 ps 1.3 lc rgb "blue"       t "m=0", \
     f1 u ($1+0.01):8:9 w yerrorbars pt 9 ps 1.4 lc rgb "dark-green" t "m=+1"

unset arrow 1; unset label 1
set output "gevp_perm_F_L2_Nf2g8_claude.png"
set title "F l=1, L=2 -- explicit per-m GEVP: m=0 block is SPURIOUS (breaks degeneracy), Nf2 g8"
set xrange [0.1:1.3]
set yrange [0.8:1.7]
set arrow 1 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2" at 1.15,sqrt(2)+0.02 tc rgb "gray30"
f2 = "permF_L2_Nf2g8_claude.dat"
plot f2 u ($1-0.01):4:5 w yerrorbars pt 7 ps 1.3 lc rgb "red"        t "m=-1 (~1.31)", \
     f2 u ($1+0.00):6:7 w yerrorbars pt 5 ps 1.4 lc rgb "blue"       t "m=0 (spurious ~1.12)", \
     f2 u ($1+0.01):8:9 w yerrorbars pt 9 ps 1.4 lc rgb "dark-green" t "m=+1 (~1.31)"
