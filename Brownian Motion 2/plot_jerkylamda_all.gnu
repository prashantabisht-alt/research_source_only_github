set multiplot layout 2,2 title "Jerk-Lambda (m x¨ + γ x˙ + λ x‴ = ξ), Gaussian noise"

# --- Panel 1: Acceleration PDF ---
set title "Acceleration PDF  P(a)"
set xlabel "a"
set ylabel "P(a)"
plot "acceleration_dist_jerk_lambda.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "#1f77b4" title "Sim"

# --- Panel 2: Velocity PDF ---
set title "Velocity PDF  P(v)"
set xlabel "v"
set ylabel "P(v)"
plot "velocity_dist_jerk_lambda.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "#ff7f0e" title "Sim"

# --- Panel 3: Position PDF ---
set title "Position PDF  P(x)"
set xlabel "x"
set ylabel "P(x)"
plot "position_dist_jerk_lambda.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "#2ca02c" title "Sim"

# --- Panel 4: Variances vs time ---
unset logscale y   # keep this one linear so growth is clear
set title "Variance vs time"
set xlabel "t"
set ylabel "Variance"
set key top left
plot \
  "jerk_lambda_stats.txt" using 1:3  with lines lw 2 lc rgb "#1f77b4" title "Var(a)", \
  ""                       using 1:5  with lines lw 2 lc rgb "#ff7f0e" title "Var(v)", \
  ""                       using 1:7  with lines lw 2 lc rgb "#2ca02c" title "Var(x)"

unset multiplot
set output
