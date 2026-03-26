reset
set term qt
set view map
unset surface
unset clabel
set contour base

# Pick contour levels (adjust as needed)
set cntrparam levels discrete 0.002,0.004,0.006,0.008,0.010,0.012

# Extract contours from simulation
set table "sim_contours.dat"
splot "joint_distribution50.txt" u 1:2:3
unset table

# Extract contours from theory
set table "theory_contours.dat"
splot "joint_theory50.txt" u 1:2:3
unset table

# Now plot both sets of contour lines
plot "sim_contours.dat"    u 1:2 w l lc rgb "black" lw 1.2 t "Simulation", \
     "theory_contours.dat" u 1:2 w l lc rgb "red"   lw 1.6 t "Theory"
