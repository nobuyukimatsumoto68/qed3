set terminal pdfcairo color size 25cm,20cm font 'Helvetica, 34'
set output 'eig_ov1_.pdf'

set xrange [-0.0:6.0]
set yrange [-2.0:2.0]
set size ratio 2.0/3.0

set lmargin at screen 0.15
set rmargin at screen 0.97
set bmargin at screen 0.15
set tmargin at screen 1.0

set key opaque
set key box vertical width 2 height 0.1 maxcols 0 spacing 1.1
set pointsize 1.4

set xlabel 'Re ~{/Symbol l}{.6\~}'
set ylabel 'Im ~{/Symbol l}{.6\~}'
set title 'generalized eigenvalues, T=4, L=4, L_t=24'

plot "geig_DW_L4_Nt24.dat" u 2:3 ti "D_W" lw 2
replot "geig_ov_L4_Nt24.dat" u 2:3 ti "Dov" ls 8 lw 2

