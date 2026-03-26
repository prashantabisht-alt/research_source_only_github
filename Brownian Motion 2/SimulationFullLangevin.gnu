set title 'Velocity Distribution from Langevin Simulation' font ',16'
set xlabel 'Velocity v'
set ylabel 'P(v)'

set key top right
set grid

# Plot simulated data (points) and theoretical Maxwell-Boltzmann (line)
plot 'velocity_distributions1.txt' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Simulation', \
     'velocity_distributions1.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Theoretical Maxwell-Boltzmann'
