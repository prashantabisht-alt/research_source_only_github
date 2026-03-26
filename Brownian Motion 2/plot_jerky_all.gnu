set multiplot layout 1,3 title "Jerky Langevin PDFs: Simulation vs Theory" font "Arial,16"

# --- Acceleration ---
set title "Acceleration PDF"
set xlabel "a"
set ylabel "P(a)"
plot "acceleration_dist.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "blue" title "Simulation", \
     "acceleration_dist.txt" using 1:3 with lines lw 2 lc rgb "red" title "Theory"

# --- Velocity ---
set title "Velocity PDF"
set xlabel "v"
set ylabel "P(v)"
plot "velocity_dist.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "blue" title "Simulation", \
     "velocity_dist.txt" using 1:3 with lines lw 2 lc rgb "red" title "Theory"

# --- Position ---
set title "Position PDF"
set xlabel "x"
set ylabel "P(x)"
plot "position_dist.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "blue" title "Simulation", \
     "position_dist.txt" using 1:3 with lines lw 2 lc rgb "red" title "Theory"

unset multiplot
