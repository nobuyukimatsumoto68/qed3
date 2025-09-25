set terminal pdfcairo color size 25cm,20cm font 'Helvetica, 34'
set output 'eig_DW_varyL_.pdf'

set xrange [-2:16]
set yrange [-6:6]
set size ratio 2.0/3.0

set lmargin at screen 0.1
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 1.0

set key opaque
set key box vertical width 2 height 0.1 maxcols 0 spacing 1.1
set pointsize 1.4

set xlabel 'Re {/Symbol l}'
set ylabel 'Im {/Symbol l}'
set title 'D_W spectrum, T=4, L_t=24'

plot "eig_DW_L4_Nt24.dat" u 2:3 ti "L=4" lw 2
replot "eig_DW_L2_Nt24.dat" u 2:3 ti "L=2" lw 2
replot "eig_DW_L1_Nt24.dat" u 2:3 ti "L=1" lw 2