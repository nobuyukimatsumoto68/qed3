set terminal pdfcairo color size 25cm,20cm font 'Helvetica, 34'
set output 'eig_DW_varyNt_.pdf'

set xrange [-0.5:17.5]
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
set title 'D_W spectrum, T=4, L=2'

plot "eig_DW_L2_Nt48.dat" u 2:3 ti " L_t=48" ls 6 lw 2
replot "eig_DW_L2_Nt24.dat" u 2:3 ti " L_t=24" ls 2 lw 2
replot "eig_DW_L2_Nt16.dat" u 2:3 ti " L_t=16" ls 12 lw 2