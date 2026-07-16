set terminal pngcairo size 1050,680 enhanced font "Helvetica,14"
set output "gevp_fit_F_L2_Nf2g8_claude.png"
set title "F l=1, L=2 -- per-m DIAG-chi2 fit + var-avg (m=0 low), Nf2 g8"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set xrange [0.1:1.3]
set yrange [0.8:1.7]
set key top right
set arrow 100 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray40"
set label 100 "{/Symbol \326}2" at graph 0.82,first sqrt(2) offset 0,0.6 tc rgb "gray40"
set object 1 rect from 0.2,1.0901898963695136 to 0.8,1.3496475505808678 fc rgb "red" fs transparent solid 0.15 noborder front
set arrow 200 from 0.2,1.2199187234751907 to 0.8,1.2199187234751907 nohead lc rgb "red" lw 2
set label 300 "m=-1: 1.220(0.130)  {/Symbol c}^2/dof=0.66" at graph 0.03,graph 0.93 tc rgb "red"
set object 2 rect from 0.2,0.9895481431186903 to 0.8,1.1621350149783574 fc rgb "blue" fs transparent solid 0.15 noborder front
set arrow 201 from 0.2,1.0758415790485238 to 0.8,1.0758415790485238 nohead lc rgb "blue" lw 2
set label 301 "m=0: 1.076(0.086)  {/Symbol c}^2/dof=0.30" at graph 0.03,graph 0.8700000000000001 tc rgb "blue"
set object 3 rect from 0.2,1.2616855829774718 to 0.8,1.4584854333760544 fc rgb "dark-green" fs transparent solid 0.15 noborder front
set arrow 202 from 0.2,1.360085508176763 to 0.8,1.360085508176763 nohead lc rgb "dark-green" lw 2
set label 302 "m=+1: 1.360(0.098)  {/Symbol c}^2/dof=0.58" at graph 0.03,graph 0.81 tc rgb "dark-green"
set object 4 rect from 0.2,1.1399492122409352 to 0.8,1.267079364186881 fc rgb "black" fs transparent solid 0.22 noborder front
set arrow 250 from 0.2,1.203514288213908 to 0.8,1.203514288213908 nohead lc rgb "black" lw 3 dt 1
set label 260 "var-avg = 1.204(0.064)" at graph 0.03,graph 0.72 tc rgb "black"
plot \
  "permF_L2_Nf2g8_claude.dat" u ($1-0.012):4:5 w yerrorbars pt 7 ps 1.2 lc rgb "red" t "m=-1", \
  "permF_L2_Nf2g8_claude.dat" u ($1+0.000):6:7 w yerrorbars pt 5 ps 1.2 lc rgb "blue" t "m=0", \
  "permF_L2_Nf2g8_claude.dat" u ($1+0.012):8:9 w yerrorbars pt 9 ps 1.2 lc rgb "dark-green" t "m=+1"
