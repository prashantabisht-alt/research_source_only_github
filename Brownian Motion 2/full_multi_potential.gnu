# --- Axes, Title, Scale ---
set logscale y
set xlabel 'Position x'
set ylabel 'P(x)'  # Correct: you're plotting P(x) on log-scale, not ln P(x)
set title 'Full Langevin Position Distributions with Different Potentials (t = 50)'
set xrange [-15:15]
set key bottom right
set grid

# --- Theoretical Distributions ---

# Harmonic: V(x) = 0.1 x^2 ⇒ P(x) ~ exp(-0.1 x^2)
# Var = 1/(2*0.1) = 5, so normalization includes sqrt(2π*5)
f_harm(x) = exp(-0.1 * x**2) / sqrt(2.0 * pi * 5.0)

# Quartic: V(x) = 0.1 x^4 ⇒ P(x) ~ exp(-0.1 x^4), approx normalization
f_quart(x) = exp(-0.1 * x**4) / (2.0 * 1.7724538509 * gamma(0.75))

# --- Plot: Points for simulation, lines for theory ---
plot \
  'pdf_harm_gauss_full.txt' using 1:2 with points pt 7 ps 0.5 lc 'blue' title 'Gaussian (Harmonic)', \
  'pdf_harm_exp_full.txt' using 1:2 with points pt 7 ps 0.5 lc 'red' title 'Exponential (Harmonic)', \
  f_harm(x) with lines lw 2 lc 'black' title 'Theoretical Harmonic', \
  'pdf_quart_gauss_full.txt' using 1:2 with points pt 7 ps 0.5 lc 'cyan' title 'Gaussian (Quartic)', \
  'pdf_quart_exp_full.txt' using 1:2 with points pt 7 ps 0.5 lc 'orange' title 'Exponential (Quartic)', \
  f_quart(x) with lines lw 2 lc 'gray' title 'Theoretical Quartic'
