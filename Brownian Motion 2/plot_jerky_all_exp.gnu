set multiplot layout 3,1 title "Exponential Noise: PDF Comparisons"

# Acceleration
set title "Acceleration PDF"
set xlabel "a"
set ylabel "P(a)"
plot "acceleration_dist_exp.txt" using 1:2 with points pt 7 lc rgb "blue" title "Simulation", \
     "acceleration_dist_exp.txt" using 1:3 with lines lw 2 lc rgb "red" title "Gaussian Theory Ref."

# Velocity
set title "Velocity PDF"
set xlabel "v"
set ylabel "P(v)"
plot "velocity_dist_exp.txt" using 1:2 with points pt 7 lc rgb "blue" title "Simulation", \
     "velocity_dist_exp.txt" using 1:3 with lines lw 2 lc rgb "red" title "Gaussian Theory Ref."

# Position
set title "Position PDF"
set xlabel "x"
set ylabel "P(x)"
plot "position_dist_exp.txt" using 1:2 with points pt 7 lc rgb "blue" title "Simulation", \
     "position_dist_exp.txt" using 1:3 with lines lw 2 lc rgb "red" title "Gaussian Theory Ref."

unset multiplot
