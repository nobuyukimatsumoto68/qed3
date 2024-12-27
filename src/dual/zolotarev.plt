set term pdfcairo enhanced
set output 'zolotarev.pdf'
set xlabel font "Liberation Sans,16"
set ylabel font "Liberation Sans,16"
set tics font "Liberation Sans,20"
set title font "Liberation Sans,22"

set xlabel 'x'
set ylabel '|r(x)|-1'

# set logscale x
set logscale y
set title '{/Symbol e}=0.01, n=21'

plot 'zolotarev.dat' using ($1):(abs($2))
