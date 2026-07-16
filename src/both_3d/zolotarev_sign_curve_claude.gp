# zolotarev_sign_curve_claude.gp
# _claude: plot the Zolotarev sign approximation for the FINALIZED split -- ACTION op D (n=31, full window
# [lambda_min, lambda_max]) vs FORCE op Df (n_f=11, NARROWED window [2*lambda_min, lambda_max]). Same lambda_max
# + same D_W, so both are functions of the physical singular value lambda. Data: zolotarev_sign_L{1,2,4}_claude.dat
# from zolotarev_sign_curve_claude.cpp. Three views per L:
#   (a) R vs lambda (log-x): action hugs sign to Delta~1e-9; force is cruder and rolls off already at 2*lambda_min,
#   (b) |1 - R| vs lambda (log-log): the two Delta plateaus (action << force),
#   (c) dR/dlambda vs lambda (log-x): the transition slope -- force peak lower AND shifted out to 2*lambda_min
#       (fewer poles + narrower window both tame it = the reversibility-clean force).
# Window edges: lambda_min (action, solid gray), 2*lambda_min (force, dashed gray), lambda_max.
# Data columns: 1=lambda, 2=R_action, 3=R_force, 4=dRdlam_action, 5=dRdlam_force.
#
# Run: gnuplot zolotarev_sign_curve_claude.gp
#   (-> zolotarev_sign_L{L}_claude.png, zolotarev_signerr_L{L}_claude.png, zolotarev_signderiv_L{L}_claude.png)

set terminal pngcairo size 1000,680 enhanced font "Helvetica,15"

# frozen windows per L (frozen_window_claude.h); arrays 1-indexed, index 4 used for L=4.
array LMIN[4]
array LMAX[4]
LMIN[1] = 0.1;   LMAX[1] = 13.0
LMIN[2] = 0.06;  LMAX[2] = 8.0
LMIN[4] = 0.008; LMAX[4] = 5.0
WINF = 2.0   # force lambda_min = WINF * lambda_min

array NACT[4]   # per-L action pole count (hasenbusch_naction): 25 (L1/L2), 31 (L4)
NACT[1] = 25; NACT[2] = 25; NACT[4] = 31

cA  = 'black'     # action (reference)
cF  = '#1f77b4'   # force n_f=11, window w=2 (2*lambda_min)
cF1 = '#d62728'   # force n_f=11, window w=1 (full, lambda_min)

do for [L in "1 2 4"] {
  Li  = L + 0
  dat = sprintf("zolotarev_sign_L%s_claude.dat", L)
  lmn = LMIN[Li]
  lmf = WINF * LMIN[Li]
  lmx = LMAX[Li]
  nact = NACT[Li]
  atitle = sprintf("action n=%d", nact)

  unset logscale x               # linear x-axis
  set xrange [0 : 8.0*lmn]       # zoom to the transition region (0 .. 8*lambda_min)
  set grid
  set key box

  # window-edge verticals + labels (shared macro via reset each panel)
  set arrow 1 from lmn, graph 0 to lmn, graph 1 nohead dt 1 lc rgb "gray55"
  set arrow 2 from lmf, graph 0 to lmf, graph 1 nohead dt 3 lc rgb "gray55"
  set arrow 3 from lmx, graph 0 to lmx, graph 1 nohead dt 1 lc rgb "gray70"

  # ---------- (a) R(lambda) ----------
  set output sprintf("zolotarev_sign_L%s_claude.png", L)
  set title sprintf("Sign approx R({/Symbol l}):  action n=%d vs force n_f=11,  L=%s  ({/Symbol l}_{min}=%.3g, force {/Symbol l}_{min}=%.3g, {/Symbol l}_{max}=%.3g)", nact, L, lmn, lmf, lmx)
  set xlabel "{/Symbol l}  (singular value of D_W)"
  set ylabel "R({/Symbol l}/{/Symbol l}_{max})"
  unset logscale y
  set yrange [0:1.12]
  set key bottom right
  set label 1 "{/Symbol l}_{min}"       at lmn, 0.05 left  offset 0.5,0 tc rgb "gray40"
  set label 2 "2{/Symbol l}_{min}(force)" at lmf, 0.11 left  offset 0.5,0 tc rgb "gray40"
  plot 1.0 w l lc rgb "gray80" dt 2 notitle, \
       dat u 1:2 w l lw 3 dt 1 lc rgb cA  title atitle, \
       dat u 1:3 w l lw 2 dt 2 lc rgb cF  title "force n_f=11 (w=2, 2{/Symbol l}_{min})", \
       dat u 1:4 w l lw 2 dt 4 lc rgb cF1 title "force n_f=11 (w=1, {/Symbol l}_{min})"
  unset label 1; unset label 2

  # ---------- (b) |1 - R(lambda)| ----------
  set output sprintf("zolotarev_signerr_L%s_claude.png", L)
  set title sprintf("Sign-approx error |1 - R({/Symbol l})|:  action n=%d vs force n_f=11,  L=%s", nact, L)
  set xlabel "{/Symbol l}  (singular value of D_W)"
  set ylabel "|1 - R({/Symbol l}/{/Symbol l}_{max})|"
  set logscale y
  set yrange [1e-11:2]
  set key top right
  set label 1 "{/Symbol l}_{min}"        at lmn, 3e-11 left offset 0.5,0 tc rgb "gray40"
  set label 2 "2{/Symbol l}_{min}(force)" at lmf, 3e-10 left offset 0.5,0 tc rgb "gray40"
  plot dat u 1:(abs(1.0-$2)) w l lw 3 dt 1 lc rgb cA  title atitle, \
       dat u 1:(abs(1.0-$3)) w l lw 2 dt 2 lc rgb cF  title "force n_f=11 (w=2, 2{/Symbol l}_{min})", \
       dat u 1:(abs(1.0-$4)) w l lw 2 dt 4 lc rgb cF1 title "force n_f=11 (w=1, {/Symbol l}_{min})"
  unset label 1; unset label 2

  # ---------- (c) dR/dlambda ----------
  set output sprintf("zolotarev_signderiv_L%s_claude.png", L)
  set title sprintf("Sign-approx slope dR/d{/Symbol l}:  action n=%d vs force n_f=11,  L=%s  (force peak lower + shifted to 2{/Symbol l}_{min})", nact, L)
  set xlabel "{/Symbol l}  (singular value of D_W)"
  set ylabel "dR/d{/Symbol l}"
  unset logscale y
  set autoscale y
  set key top right
  set label 1 "{/Symbol l}_{min}"        at lmn, graph 0.05 left offset 0.5,0 tc rgb "gray40"
  set label 2 "2{/Symbol l}_{min}(force)" at lmf, graph 0.11 left offset 0.5,0 tc rgb "gray40"
  plot 0 w l lc rgb "gray80" dt 2 notitle, \
       dat u 1:5 w l lw 3 dt 1 lc rgb cA  title atitle, \
       dat u 1:6 w l lw 2 dt 2 lc rgb cF  title "force n_f=11 (w=2, 2{/Symbol l}_{min})", \
       dat u 1:7 w l lw 2 dt 4 lc rgb cF1 title "force n_f=11 (w=1, {/Symbol l}_{min})"
  unset label 1; unset label 2

  unset arrow 1; unset arrow 2; unset arrow 3
}
