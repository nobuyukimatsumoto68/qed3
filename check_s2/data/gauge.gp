set terminal pdfcairo color size 20cm,20cm font 'Helvetica, 34'
set output 'delta_s2_.pdf'

set xrange [-1.0:1.0]
set yrange [-20:100.0]
set size ratio 1.0

set ytics
set xtics

set lmargin at screen 0.15
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.9

set key opaque
set key box vertical width 2 height 0.4 spacing 1.5

set xlabel 'z'
set title 'J^t-J^t correlator on S^2'

plot "gauge_n40.dat" u 1:2 ti 'n_{max}=40' w l ls 1 lw 7
replot "gauge_n20.dat" u 1:2 ti 'n_{max}=20' w l ls 2 lw 7 dashtype 4
