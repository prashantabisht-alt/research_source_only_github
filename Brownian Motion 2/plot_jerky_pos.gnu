set xlabel "x"
set ylabel "P(x)"
plot "position_dist.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "blue" title "Simulation", \
     "position_dist.txt" using 1:3 with lines lw 2 lc rgb "red" title "Theory"

