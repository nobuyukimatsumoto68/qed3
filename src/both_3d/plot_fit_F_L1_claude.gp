set terminal pngcairo size 1050,680 enhanced font "Helvetica,14"
set output "gevp_fit_F_L1_Nf2g8_claude.png"
set title "F l=1, L=1 -- per-m DIAG-chi2 fit + var-avg, Nf2 g8"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:1.3]
set yrange [1.15:1.5]
set key top right
set arrow 100 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray40"
set label 100 "{/Symbol \326}2" at graph 0.82,first sqrt(2) offset 0,0.6 tc rgb "gray40"
set object 1 rect from 0.2,1.3434145838982123 to 0.8,1.4374569965665438 fc rgb "red" fs transparent solid 0.15 noborder front
set arrow 200 from 0.2,1.390435790232378 to 0.8,1.390435790232378 nohead lc rgb "red" lw 2
set label 300 "m=-1: 1.390(0.047)  {/Symbol c}^2/dof=0.04" at graph 0.03,graph 0.93 tc rgb "red"
set object 2 rect from 0.2,1.3210483399111383 to 0.8,1.4116419385861485 fc rgb "blue" fs transparent solid 0.15 noborder front
set arrow 201 from 0.2,1.3663451392486434 to 0.8,1.3663451392486434 nohead lc rgb "blue" lw 2
set label 301 "m=0: 1.366(0.045)  {/Symbol c}^2/dof=0.00" at graph 0.03,graph 0.8700000000000001 tc rgb "blue"
set object 3 rect from 0.2,1.2248362670131274 to 0.8,1.3293758202490555 fc rgb "dark-green" fs transparent solid 0.15 noborder front
set arrow 202 from 0.2,1.2771060436310915 to 0.8,1.2771060436310915 nohead lc rgb "dark-green" lw 2
set label 302 "m=+1: 1.277(0.052)  {/Symbol c}^2/dof=0.08" at graph 0.03,graph 0.81 tc rgb "dark-green"
set object 4 rect from 0.2,1.3216491369803607 to 0.8,1.37769907272406 fc rgb "black" fs transparent solid 0.22 noborder front
set arrow 250 from 0.2,1.3496741048522103 to 0.8,1.3496741048522103 nohead lc rgb "black" lw 3 dt 1
set label 260 "var-avg = 1.350(0.028)" at graph 0.03,graph 0.72 tc rgb "black"
plot \
  "permF_L1_Nf2g8_claude.dat" u ($1-0.012):4:5 w yerrorbars pt 7 ps 1.2 lc rgb "red" t "m=-1", \
  "permF_L1_Nf2g8_claude.dat" u ($1+0.000):6:7 w yerrorbars pt 5 ps 1.2 lc rgb "blue" t "m=0", \
  "permF_L1_Nf2g8_claude.dat" u ($1+0.012):8:9 w yerrorbars pt 9 ps 1.2 lc rgb "dark-green" t "m=+1"
