set multiplot layout 2,1 title "Velocity and Position Distributions from Full Langevin Simulation"

# --- First plot: Velocity distribution ---
set xlabel 'Velocity v'
set ylabel 'P(v)'
set title 'Velocity Distribution vs Maxwell-Boltzmann Theory'
set grid
plot 'velocity_distributions20.txt' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Simulation', \
     'velocity_distributions20.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Maxwell-Boltzmann Theory'

# --- Second plot: Position distribution ---
set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Position Distribution vs Diffusive Gaussian Theory'
set grid
plot 'position_distributions20.txt' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Simulation', \
     'position_distributions20.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Gaussian Theory'

unset multiplot
