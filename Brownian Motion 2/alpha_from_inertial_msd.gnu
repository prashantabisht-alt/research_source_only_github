set logscale x
set xlabel "t / τI"
set ylabel "α(t) = d ln(MSD)/d ln t"

plot "alpha_active_taupratio_1e-3two.txt" u 1:2 w l lw 2 t "τp = τI/1000", \
     "alpha_active_taupratio_1e3two.txt"  u 1:2 w l lw 2 t "τp = 1000 τI", \
     6 w l lw 1 dt 2 lc rgb "black" t "y=6", \
     5 w l lw 1 dt 2 lc rgb "black" t "y=5", \
     4 w l lw 1 dt 2 lc rgb "black" t "y=4", \
     3 w l lw 1 dt 2 lc rgb "black" t "y=3"
