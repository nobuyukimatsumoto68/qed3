#!/bin/bash
# plot_WI_wgp_claude.sh
# Usage: ./plot_WI_wgp_claude.sh [logfile]
# Extracts running-average Ward identity data from the check log and plots with gnuplot.

LOG=${1:-current_check.log}
DAT=WI_running_avg_claude.dat

grep "# hit" "$LOG" \
  | awk '{split($4,a,"="); print $3, a[2]}' \
  > "$DAT"

echo "# Extracted $(wc -l < "$DAT") points from $LOG -> $DAT"

gnuplot -p <<'EOF'
set xlabel "hit"
set ylabel "|running avg of div(w)|"
set title "Stochastic Ward identity running average"
set logscale y
set grid
plot "WI_running_avg_claude.dat" using 1:(abs($2)) with lines title "|running avg|"
EOF
