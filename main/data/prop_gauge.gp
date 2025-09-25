set terminal pdfcairo color size 20cm,20cm font 'Helvetica, 34'
set output 'jtjt_gaussian_.pdf'

set xrange [0.:12.0]
set size square

set logscale y

set lmargin at screen 0.15
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.9

set yrange [0.00001:1000]
set format y "10^{%T}"
# set ytic(12345)

set xlabel 't'
set title 'J^t-J^t correlator with 1/g^2=20, T=12, L_t=120'

plot "gauge_exact.dat" u 1:2 ti "exact" w l ls 3 lw 3 lt rgb "black" dashtype 7

set key opaque
set key box vertical width 1.2 height 0.5 maxcols 0 spacing 1.1
set pointsize 2

replot "prop_gauge_gsq0.050000_L8_Nt120_T12.000000.dat" u (12.*$0/120):(20.*$1) ti "L=8" w p ls 6 lw 3
replot "prop_gauge_gsq0.050000_L4_Nt120_T12.000000.dat" u (12.*$0/120):(20.*$1) ti "L=4" w p ls 4 lw 3
replot "prop_gauge_gsq0.050000_L2_Nt120_T12.000000.dat" u (12.*$0/120):(20.*$1) ti "L=2" w p ls 2 lw 3
replot "prop_gauge_gsq0.050000_L1_Nt120_T12.000000.dat" u (12.*$0/120):(20.*$1) ti "L=1" w p ls 1 lw 3
