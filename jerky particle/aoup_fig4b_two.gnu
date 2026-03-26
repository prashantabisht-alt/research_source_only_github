set term qt size 1000,650
set logscale x
set xlabel "s = t/τ_I"
set ylabel "α(s) = d ln(MSD)/d ln t"
set grid
set key top right spacing 1.2

plot \
  "alpha_taupratio_1e3two.txt"                u 1:2 w l lw 4 dt 2 lc rgb "red"   t "τ_p/τ_I=1e3 (sim)", \
  "alpha_theory_taupratio_1e3_overlay_numeric.txt"  u 1:2 w l lw 3 lc rgb "black" t "τ_p/τ_I=1e3 (theory)", \
  6 w l dt 3 lw 1.5 lc rgb "black"  t "α=6", \
  5 w l dt 3 lw 1.5 lc rgb "black"  t "α=5", \
  4 w l dt 3 lw 1.5 lc rgb "purple" t "α=4", \
  3 w l dt 3 lw 1.5 lc rgb "green"  t "α=3"
