set term qt size 1000,650
set logscale x
set xlabel "s = t/τ_I"
set ylabel "α(s) = d ln(MSD)/d ln t"
set grid
plot "alpha_taupratio_1e-3two.txt" u 1:2 w l lw 2 lc rgb "red"  t "τ_p/τ_I=1e-3: 6→5→3", \
     "alpha_taupratio_1e3two.txt"  u 1:2 w l lw 2 lc rgb "blue" t "τ_p/τ_I=1e3 : 6→4→3", \
     6 dt 3 lc rgb "black"  t "α=6", \
     5 dt 3 lc rgb "black" t "α=5", \
     4 dt 3 lc rgb "purple" t "α=4", \
     3 dt 3 lc rgb "green"  t "α=3"
