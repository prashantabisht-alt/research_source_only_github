set logscale x
set xlabel "t/τ_p"
set ylabel "α(t) = d ln MSD / d ln t"
set yrange [4.5:6.5]
plot "alpha_persistent.txt" u 1:2 w l lw 2 t "simulation α(t)", \
     6.0 lc rgb "gray" dt 2 t "t^6", \
     5.0 lc rgb "gray" dt 2 t "t^5"
