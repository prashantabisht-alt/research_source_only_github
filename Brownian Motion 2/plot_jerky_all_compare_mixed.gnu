set multiplot layout 1,3 title "Gaussian vs Exponential Noise (log scale)" font ",16"

# --- Acceleration ---
set xlabel "a"
set ylabel "P(a)"
plot "acceleration_dist_mixed_gaussian.txt"    u 1:2 w l lw 2 lc rgb "blue" title "Gaussian", \
     "acceleration_dist_mixed_exponential.txt" u 1:2 w l lw 2 lc rgb "red"  title "Exponential"

# --- Velocity ---
set xlabel "v"
set ylabel "P(v)"
plot "velocity_dist_mixed_gaussian.txt"    u 1:2 w l lw 2 lc rgb "blue" title "Gaussian", \
     "velocity_dist_mixed_exponential.txt" u 1:2 w l lw 2 lc rgb "red"  title "Exponential"

# --- Position ---
set xlabel "x"
set ylabel "P(x)"
plot "position_dist_mixed_gaussian.txt"    u 1:2 w l lw 2 lc rgb "blue" title "Gaussian", \
     "position_dist_mixed_exponential.txt" u 1:2 w l lw 2 lc rgb "red"  title "Exponential"

unset multiplot
