set term pdfcairo enhanced
set output 'stepsize_scaling.pdf'
# set term pdf enhanced
# set output 'stepsize_scaling.svg'
set xlabel font "Liberation Sans,16"
set ylabel font "Liberation Sans,16"
set tics font "Liberation Sans,20"
set title font "Liberation Sans,22"

set xlabel '{/Symbol t}'
set ylabel '{/Symbol D} H'

set logscale x
set logscale y

plot 'scaling.dat' using (1.0/$1):(abs($2))
replot 100.0 * x**2
