set term qt size 900,550
set xlabel "time t"
set ylabel "position x(t)"
set grid
plot "Brownian_single.txt" u 1:2 w l lw 1.6 lc rgb "blue" t "1D Brownian trajectory"
