set terminal pngcairo size 1050,680 enhanced font "Helvetica,14"
set output "gevp_fit_F2_L1_Nf2g8_claude.png"
set title "F^2 0++, L=1 -- DIAG-chi2 fit, Nf2 g8"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:1.3]
set yrange [-0.2:3.5]
set key top right
set arrow 100 from 0.1,2*sqrt(2) to 1.3,2*sqrt(2) nohead dt 2 lc rgb "gray40"
set label 100 "2{/Symbol \326}2" at graph 0.82,first 2*sqrt(2) offset 0,0.6 tc rgb "gray40"
set object 1 rect from 0.2,2.7008828760161463 to 0.6,2.994061332702763 fc rgb "blue" fs transparent solid 0.15 noborder front
set arrow 200 from 0.2,2.8474721043594546 to 0.6,2.8474721043594546 nohead lc rgb "blue" lw 2
set label 300 "0++: 2.847(0.147)  {/Symbol c}^2/dof=0.93" at graph 0.03,graph 0.93 tc rgb "blue"
plot \
  "permF2_L1_eff_Nf2g8_claude.dat" u ($1-0.012):4:5 w yerrorbars pt 7 ps 1.2 lc rgb "blue" t "0++"
