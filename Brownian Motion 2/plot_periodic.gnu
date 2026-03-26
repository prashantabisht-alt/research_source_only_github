set logscale y
set xlabel 'Position x'
set ylabel 'P(x)'  # You're plotting P(x), not log(P)
set title 'Overdamped Position Distributions in Periodic Potential (t = 500)'
set xrange [-15:15]
set key bottom right
set grid

# Theoretical Boltzmann-like distribution in periodic potential
f(x) = exp(-2.0 * cos(2.0 * x)) / (2.0 * pi)  # Approximate normalization

plot \
  'pos_dist_periodic_gaussian.txt' using 1:2 with points pt 7 ps 0.5 lc 'blue' title 'Gaussian Noise', \
  'pos_dist_periodic_exponential.txt' using 1:2 with points pt 7 ps 0.5 lc 'red' title 'Exponential Noise', \
  f(x) with lines lw 2 lc 'black' title 'Theoretical Periodic'
