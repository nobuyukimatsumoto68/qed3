set term pdfcairo enhanced
set output 'cost_scaling.pdf'
set xlabel font "Liberation Sans,16"
set ylabel font "Liberation Sans,16"
set tics font "Liberation Sans,20"
set title font "Liberation Sans,22"

set xlabel 'volume'
set ylabel 'sec'
set title "10 configs"

set logscale x
set logscale y

plot 'cost_scaling.dat' using ($1):(abs($2))
replot 0.0001 * x**2
