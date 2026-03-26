set title "Stationary Velocity Distribution: Gaussian vs. Exponential Noise" font ",16"
set xlabel "Velocity (v)"
set ylabel "Probability Density P(v)"
set key top right
set grid

# --- Plot Command ---
# The theoretical curve is read from the third column of the Gaussian data file.
# We plot both simulated distributions to see if they match the theory.

plot 'dist_gaussian.txt' using 1:3 with lines lw 2 lc 'black' title 'Theoretical Maxwell-Boltzmann', \
     'dist_gaussian.txt' using 1:2 with points pt 7 ps 0.5 lc 'blue' title 'Simulation (Gaussian Noise)', \
     'dist_exponential.txt' using 1:2 with points pt 5 ps 0.5 lc 'red' title 'Simulation (Exponential Noise)'
