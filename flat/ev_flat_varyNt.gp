set terminal pdfcairo color size 25cm,20cm font 'Helvetica, 34'
set output 'flat_spectrum_varyNt_.pdf'

set xrange [-0.5:17.5]
set yrange [-6:6]
set size ratio 2.0/3.0

set lmargin at screen 0.1
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 1.0

set key opaque
set key box vertical width 1 height 0.1 maxcols 0 spacing 1.1
set pointsize 1.4

set xlabel 'Re {/Symbol l}_{flat}'
set ylabel 'Im {/Symbol l}_{flat}'
set title 'D_{W} spectrum (flat limit from T=4, L=2)'

plot "flat_L2_Lt48_1.dat" u 1:2 ti "L_t=48" ls 6 lw 2
replot "flat_L2_Lt48_2.dat" u 1:2 notitle ls 6 lw 2
replot "flat_L2_Lt24_1.dat" u 1:2 ti "L_t=24" ls 2 lw 2
replot "flat_L2_Lt24_2.dat" u 1:2 notitle ls 2 lw 2
replot "flat_L2_Lt16_1.dat" u 1:2 ti "L_t=16" ls 12 pt 3 lw 2
replot "flat_L2_Lt16_2.dat" u 1:2 notitle ls 12 pt 3 lw 2

