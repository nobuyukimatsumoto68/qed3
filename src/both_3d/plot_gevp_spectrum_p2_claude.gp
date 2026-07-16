# plot_gevp_spectrum_p2_claude.gp
# GEVP spectra (per-state effmass vs t) for the phase-2 setup, Nf2 g8. Three panels/outputs:
#  (1) F l=1, L=1, TRIANGLE-only: the 3 m-states (should be ~degenerate = T_1 irrep).
#  (2) F l=1, L=2, FULL basis: the l=1 tower -- tests whether the low Delta~1.12 mode is a real
#      plateau or a spurious multi-shape GEVP artifact.
#  (3) F^2 0++, L=1, FULL basis: vacuum(Delta~0) + physical 0++(2sqrt2) + excited.
# cols: t, ground(2,3), then s0,s1,... at (4,5),(6,7),...  (s_last = ground = lightest).
set terminal pngcairo size 1000,660 enhanced font "Helvetica,14"
set xlabel "t"
set ylabel "{/Symbol D}_{eff}"
set grid
set key top right

# ---- (1) F l=1 L1 triangle: 3 m-states ----
set output "gevp_spectrum_p2_F_L1tri_Nf2g8_claude.png"
set title "F l=1 GEVP spectrum -- L=1, triangle-only, Nf2 g8 (3 m-states)"
set xrange [0.1:1.5]
set yrange [1.0:1.7]
set arrow 1 from 0.1,sqrt(2) to 1.5,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2 (free, continuum)" at 0.9,sqrt(2)+0.02 tc rgb "gray30"
f1 = "spec_F_L1tri_Nf2g8_claude.dat"
plot f1 u ($1-0.01):4:5 w yerrorbars pt 7 ps 1.2 lc rgb "red"        t "m-state a", \
     f1 u ($1+0.00):6:7 w yerrorbars pt 5 ps 1.2 lc rgb "blue"       t "m-state b", \
     f1 u ($1+0.01):8:9 w yerrorbars pt 9 ps 1.3 lc rgb "dark-green" t "m-state c"

# ---- (2) F l=1 L2 full: tower ----
unset arrow 1; unset label 1
set output "gevp_spectrum_p2_F_L2full_Nf2g8_claude.png"
set title "F l=1 GEVP spectrum -- L=2, full 5-shape, Nf2 g8 (is Delta~1.1 real?)"
set xrange [0.1:1.3]
set yrange [-1.0:5.5]
set arrow 1 from 0.1,sqrt(2) to 1.3,sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "{/Symbol \326}2" at 1.1,sqrt(2)+0.15 tc rgb "gray30"
f2 = "spec_F_L2full_Nf2g8_claude.dat"
plot f2 u ($1+0.00):14:15 w yerrorbars pt 7 ps 1.3 lc rgb "red"        t "ground (lightest)", \
     f2 u ($1+0.00):12:13 w yerrorbars pt 5 ps 1.1 lc rgb "blue"       t "state 1", \
     f2 u ($1+0.00):10:11 w yerrorbars pt 9 ps 1.1 lc rgb "dark-green" t "state 2", \
     f2 u ($1+0.00):8:9   w yerrorbars pt 11 ps 1.1 lc rgb "orange"    t "state 3", \
     f2 u ($1+0.00):6:7   w yerrorbars pt 13 ps 1.1 lc rgb "purple"    t "state 4", \
     f2 u ($1+0.00):4:5   w yerrorbars pt 6 ps 1.1 lc rgb "black"      t "state 5"

# ---- (3) F^2 0++ L1 full: tower ----
unset arrow 1; unset label 1
set output "gevp_spectrum_p2_F2_L1full_Nf2g8_claude.png"
set title "F^2 0++ GEVP spectrum -- L=1, full 5-shape, Nf2 g8"
set xrange [0.1:1.5]
set yrange [-0.5:5.5]
set arrow 1 from 0.1,2*sqrt(2) to 1.5,2*sqrt(2) nohead dt 2 lc rgb "gray30"
set label 1 "2{/Symbol \326}2 (two-photon)" at 0.9,2*sqrt(2)+0.15 tc rgb "gray30"
set arrow 2 from 0.1,0 to 1.5,0 nohead dt 3 lc rgb "dark-green"
f3 = "spec_F2_L1full_Nf2g8_claude.dat"
plot f3 u ($1+0.00):10:11 w yerrorbars pt 7 ps 1.3 lc rgb "dark-red"   t "ground (vacuum mode)", \
     f3 u ($1+0.00):8:9   w yerrorbars pt 5 ps 1.2 lc rgb "blue"       t "physical 0++", \
     f3 u ($1+0.00):6:7   w yerrorbars pt 9 ps 1.1 lc rgb "orange"     t "excited", \
     f3 u ($1+0.00):4:5   w yerrorbars pt 6 ps 1.1 lc rgb "black"      t "state 3"
