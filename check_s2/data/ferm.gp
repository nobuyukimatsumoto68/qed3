set terminal pdfcairo color size 20cm,20cm font 'Helvetica, 34'
set output 'prop_s2_.pdf'

set xrange [0.001:pi]
set yrange [0.01:40.0]
set size ratio 1.0

set ytics
set xtics

set logscale xy

set lmargin at screen 0.15
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.9

set key opaque
set key box vertical width 2 height 0.1 maxcols 0 spacing 1.1

set xlabel '{/Symbol q}'
set title 'Dirac propagator on S^2'

plot "ferm_n200.dat" u 1:2 ti 'n_{max}=200' w l ls 1 lw 5
replot "ferm_n50.dat" u 1:2 ti 'n_{max}=50' w l ls 2 lw 5 dashtype 4
replot 1.0/sin(x/2.0)/(4.0*pi) ti 'exact' ls 3 lw 3 lt rgb "black" dashtype 7
