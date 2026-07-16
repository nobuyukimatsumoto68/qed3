# plot_tseries_1pt_claude.gp
# Monte-Carlo time series of the F l=1 one-point function magnitude (norm over the 12 l=1 ops)
# vs config index, for the two suspicious ensembles (Nf4/Nf6 g2 L1). Stationary series, no long
# frozen stretch (longest = 5 identical configs = normal HMC rejections at ~10-13%).
#   gnuplot plot_tseries_1pt_claude.gp -> gevp_tseries_1pt_claude.png
# cols of tseries_1pt_<tag>_claude.dat: 1=k 2=F_tri_l0(~0) 3=F_tri_l1m0 4=norm_l1(12ops)

set terminal pngcairo size 1200,640 enhanced font "Helvetica,14"
set output "gevp_tseries_1pt_claude.png"
set title "F l=1 one-point-function magnitude vs MC time (g2, L1) -- stationary, no stuck config"
set xlabel "config index k"
set ylabel "|F_{l=1}| (norm of 12 l=1 ops)"
set grid
set key top right
set yrange [0:0.25]

plot \
  "tseries_1pt_Nf4_g2_claude.dat" u 1:4 w points pt 7 ps 0.4 lc rgb "blue"       t "Nf4 g2 (5899 cfg)", \
  "tseries_1pt_Nf6_g2_claude.dat" u 1:4 w points pt 5 ps 0.4 lc rgb "dark-green" t "Nf6 g2 (5538 cfg)"
