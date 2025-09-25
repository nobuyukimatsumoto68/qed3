# set terminal postscript color enhanced size 20cm,20cm # font "cmr,10" fontscale 0.7 size 3.5in,2.4in
set terminal pdfcairo color size 20cm,20cm font 'Helvetica, 34'
set output 'ferm_prop_w_ov_.pdf'

set xrange [0.:12.]
# set yrange [-1:20.0]
set size ratio 1.0

set logscale y
set yrange [0.0001:20.0]

set lmargin at screen 0.15
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.9

set xlabel 't'
set title 'D_{lat} propagators with L=1, T=12, L_t=168'

set key opaque
set key box vertical width 2.0 height 0.5 maxcols 0 spacing 1.1
set pointsize 2

plot "../prop_temporal_L1_Nt168_T12.000000.dat" u 1:($2) every 2:2 ti 'D_W' w p ls 6 pt 3 lw 3
replot "../ov_prop_temporal_L1_Nt168_T12.000000.dat" u 1:($2) every 2:2 ti 'D_{ov}' w p ls 1 lw 3
