set multiplot layout 3,1 title "Gaussian vs Exponential Noise: PDFs"

# Acceleration comparison
set title "Acceleration PDF"
set xlabel "a"
set ylabel "P(a)"
plot "acceleration_dist.txt" using 1:2 with points pt 7 lc rgb "blue" title "Gaussian Sim", \
     "acceleration_dist_exp.txt" using 1:2 with points pt 7 lc rgb "green" title "Exponential Sim", \
     "acceleration_dist.txt" using 1:3 with lines lw 2 lc rgb "red" title "Theory"

# Velocity comparison
set title "Velocity PDF"
set xlabel "v"
set ylabel "P(v)"
plot "velocity_dist.txt" using 1:2 with points pt 7 lc rgb "blue" title "Gaussian Sim", \
     "velocity_dist_exp.txt" using 1:2 with points pt 7 lc rgb "green" title "Exponential Sim", \
     "velocity_dist.txt" using 1:3 with lines lw 2 lc rgb "red" title "Theory"

# Position comparison
set title "Position PDF"
set xlabel "x"
set ylabel "P(x)"
plot "position_dist.txt" using 1:2 with points pt 7 lc rgb "blue" title "Gaussian Sim", \
     "position_dist_exp.txt" using 1:2 with points pt 7 lc rgb "green" title "Exponential Sim", \
     "position_dist.txt" using 1:3 with lines lw 2 lc rgb "red" title "Theory"

unset multiplot
