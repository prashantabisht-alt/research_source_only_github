# ====================
# 1. Velocity Distribution
# ====================
set xlabel 'Velocity v'
set ylabel 'P(v)'
set title 'Velocity Distribution (Maxwell–Boltzmann)'
set key top right
plot 'velocity_distribution_potential.txt' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Simulation', \
     'velocity_distribution_potential.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Maxwell–Boltzmann Theory'

# ====================
# 2. Position Distribution
# ====================
set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Position Distribution (Boltzmann in Harmonic Potential)'
set key top right
plot 'position_distribution_potential1.txt' using 1:2 with points pt 7 ps 0.5 lc rgb 'green' title 'Simulation', \
     'position_distribution_potential1.txt' using 1:3 with lines lw 2 lc rgb 'orange' title 'Boltzmann Theory'

unset multiplot
