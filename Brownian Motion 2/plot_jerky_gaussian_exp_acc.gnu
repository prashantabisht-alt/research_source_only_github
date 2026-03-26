set title "Acceleration PDF"
set xlabel "a"
set ylabel "P(a)"
plot "acceleration_dist.txt" using 1:2 with points pt 3 lc rgb "blue" title "Gaussian Sim", \
     "acceleration_dist_exp.txt" using 1:2 with points pt 3 lc rgb "green" title "Exponentia>
     "acceleration_dist.txt" using 1:3 with lines lw 2 lc rgb "red" title "Theory"
