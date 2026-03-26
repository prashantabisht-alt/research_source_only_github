set term qt size 900,600
set logscale x
set xlabel "s = t/τp"
set ylabel "α(s)"
set grid
plot "alpha_persistent.txt" u 1:2 w l lw 2 t "simulation α(s)", \
     6 lc rgb "gray" dt 2 t "α=6", \
     5 lc rgb "black" dt 2 t "α=5"
