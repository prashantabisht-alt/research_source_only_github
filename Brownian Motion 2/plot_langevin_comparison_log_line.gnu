set logscale y
set title "Stationary Velocity Distribution: Gaussian vs. Exponential Noise" font ",16"
set xlabel "Velocity (v)"
set ylabel "Probability Density P(v)"
set key top right
set grid

plot \
  'dist_gaussian.txt' using 1:3 with lines lw 2 lc 'black' title 'Theoretical Maxwell-Boltzmann', \
  'dist_gaussian.txt' using 1:2 with linespoints pt 7 ps 0.5 lw 2 lc 'blue' title 'Simulation (Gaussian Noise)', \
  'dist_exponential.txt' using 1:2 with linespoints pt 5 ps 0.5 lw 2 lc 'red' title 'Simulation (Exponential Noise)'

