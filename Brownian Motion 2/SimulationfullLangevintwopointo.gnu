set multiplot layout 2,1 title "Full Langevin Simulation: Velocity & Position Distributions"

# ---------- Top plot: Velocity ----------
set title "Velocity Distribution (Maxwell-Boltzmann check)"
set xlabel "Velocity v"
set ylabel "P(v)"
set grid
set key top right
plot 'velocity_distributions20.txt' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Simulation', \
     'velocity_distributions20.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Theoretical Maxwell-Boltzmann'

# ---------- Bottom plot: Position ----------
set title "Position Distribution (Diffusion check)"
set xlabel "Position x"
set ylabel "P(x)"
set grid
set key top right
plot 'position_distributions20.txt' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Simulation', \
     'position_distributions20.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Theoretical Diffusive Gaussian'

unset multiplot
