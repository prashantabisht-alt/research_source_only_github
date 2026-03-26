set term qt size 900,550
set xlabel "time t"
set ylabel "position x(t)"
set grid
plot "Brownian_single.txt" every ::0::50 u 1:2 w lp lw 1.6 pt 7 ps 1.5 lc rgb "blue" t "first 50 steps"

