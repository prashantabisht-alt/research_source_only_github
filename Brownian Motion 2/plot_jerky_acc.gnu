set title "Acceleration PDF"
set xlabel "a"
set ylabel "P(a)"
plot "acceleration_dist.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "blue" title "Simulation", \
     "acceleration_dist.txt" using 1:3 with lines lw 2 lc rgb "red" title "Theory"
