# wilson_lowmode_claude.gp
# _claude: plot the Wilson low-mode (lambda_min = smallest singular value of D_W) distribution and
# time series for L=1,2,4 from eig_wilson_lowmode_claude.cu (gsq=8, Nf=2).
# Page 1 = 3 SEPARATE histogram panels (one per L, each with its own binning/range, since the L's live
# on very different scales). Page 2 = lambda_min history, overlaid.
# Series style (color-blind: distinct COLOR + MARKER/PATTERN): L=1 red circle (lc 1 pt 7),
# L=2 blue square (lc 3 pt 5), L=4 dark-green triangle (lc rgb '#009900' pt 9).
# Run:  gnuplot wilson_lowmode_claude.gp   ->  wilson_lowmode_claude.pdf

# Nf / gsq / terminal selectable:  gnuplot wilson_lowmode_claude.gp   (Nf=2, gsq=8, PDF)
#   or  gnuplot -e "NF=6; GSQ=12" ...            (other ensemble, PDF)
#   or  gnuplot -e "TERM='png'" ...              (PNG per page, for markdown embedding)
if (!exists("NF"))   NF   = 2
if (!exists("GSQ"))  GSQ  = 8
if (!exists("TERM")) TERM = 'pdf'

if (TERM eq 'png') {
  set terminal pngcairo color size 1400,1700 font 'Helvetica,16'
  # PNG: output is set PER PAGE below (_lmin / _lmax / _hist).
} else {
  set terminal pdfcairo color size 25cm,28cm font 'Helvetica,18'
  set output sprintf('wilson_lowmode_Nf%d_gsq%g_claude.pdf', NF, GSQ)
}

# NOTE: the driver formats gsq via std::to_string(double) -> "8.000000" / "12.000000" in the filename.
f1 = sprintf('wilson_lowmode_Nf%d_gsq%.6f_L1_claude.dat', NF, GSQ)
f2 = sprintf('wilson_lowmode_Nf%d_gsq%.6f_L2_claude.dat', NF, GSQ)
f4 = sprintf('wilson_lowmode_Nf%d_gsq%.6f_L4_claude.dat', NF, GSQ)

# --- per-L normalization + per-L bin width (each panel binned on its OWN range) ---
nbins = 30
bin(x,w) = w*(floor(x/w)+0.5)

stats f1 u 2 nooutput ; N1 = STATS_records ; w1 = (STATS_max-STATS_min)/nbins
stats f2 u 2 nooutput ; N2 = STATS_records ; w2 = (STATS_max-STATS_min)/nbins
stats f4 u 2 nooutput ; N4 = STATS_records ; w4 = (STATS_max-STATS_min)/nbins

stats f1 u 3 nooutput ; w1m = (STATS_max-STATS_min)/nbins
stats f2 u 3 nooutput ; w2m = (STATS_max-STATS_min)/nbins
stats f4 u 3 nooutput ; w4m = (STATS_max-STATS_min)/nbins

# ================= page 1: lambda_min distribution, one panel per L =================
if (TERM eq 'png') { set output sprintf('wilson_lowmode_Nf%d_gsq%g_lmin.png', NF, GSQ) }
set style fill pattern border -1
set grid
set ylabel 'fraction of configs'

set multiplot layout 3,1 title sprintf('Wilson low-mode distribution, gsq=%g, Nf=%d', GSQ, NF) font 'Helvetica,22'

  set xlabel '{/Symbol l}_{min}(D_W)   L=1'
  set boxwidth w1
  plot f1 u (bin($2,w1)):(1.0/N1) smooth freq with boxes lc 1             fs pattern 1 title 'L=1'

  set xlabel '{/Symbol l}_{min}(D_W)   L=2'
  set boxwidth w2
  plot f2 u (bin($2,w2)):(1.0/N2) smooth freq with boxes lc 3             fs pattern 2 title 'L=2'

  set xlabel '{/Symbol l}_{min}(D_W)   L=4'
  set boxwidth w4
  plot f4 u (bin($2,w4)):(1.0/N4) smooth freq with boxes lc rgb '#009900' fs pattern 5 title 'L=4'

unset multiplot

# ================= page 2: lambda_MAX distribution, one panel per L =================
if (TERM eq 'png') { set output sprintf('wilson_lowmode_Nf%d_gsq%g_lmax.png', NF, GSQ) }
set multiplot layout 3,1 title sprintf('Wilson {/Symbol l}_{max} distribution, gsq=%g, Nf=%d', GSQ, NF) font 'Helvetica,22'

  set xlabel '{/Symbol l}_{max}(D_W)   L=1'
  set boxwidth w1m
  plot f1 u (bin($3,w1m)):(1.0/N1) smooth freq with boxes lc 1             fs pattern 1 title 'L=1'

  set xlabel '{/Symbol l}_{max}(D_W)   L=2'
  set boxwidth w2m
  plot f2 u (bin($3,w2m)):(1.0/N2) smooth freq with boxes lc 3             fs pattern 2 title 'L=2'

  set xlabel '{/Symbol l}_{max}(D_W)   L=4'
  set boxwidth w4m
  plot f4 u (bin($3,w4m)):(1.0/N4) smooth freq with boxes lc rgb '#009900' fs pattern 5 title 'L=4'

unset multiplot

# ================= page 3: lambda_min history (overlaid) =================
if (TERM eq 'png') { set output sprintf('wilson_lowmode_Nf%d_gsq%g_hist.png', NF, GSQ) }
set xlabel 'trajectory k'
set ylabel '{/Symbol l}_{min}(D_W)'
set title sprintf('Wilson low-mode history, gsq=%g, Nf=%d', GSQ, NF)
set key top right box opaque
set logscale y
set pointsize 1.0

plot \
  f1 u 1:2 with points lc 1             pt 7 title 'L=1', \
  f2 u 1:2 with points lc 3             pt 5 title 'L=2', \
  f4 u 1:2 with points lc rgb '#009900' pt 9 title 'L=4'

unset logscale y
