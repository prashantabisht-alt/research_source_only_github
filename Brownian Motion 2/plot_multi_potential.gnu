# --- Log Scale for y-axis to highlight tails ---
set logscale y

# --- Axis Labels and Title ---
set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Overdamped Position Distributions with Different Potentials (t = 10)'

# --- Axis Range and Legend ---
set xrange [-10:10]
set key bottom right
set grid

# --- Theoretical Distributions ---

# Harmonic: P(x) ~ exp(-0.5 * x^2)
f_harm(x) = exp(-0.5 * x**2) / sqrt(2.0 * pi)

# Quartic: P(x) ~ exp(-0.1 * x^4), approx normalized
f_quart(x) = exp(-0.1 * x**4) / (2.0 * 1.7724538509 * gamma(0.75))

# --- Plot Command ---
plot \
  'pdf_gauss_harmonic.txt' using 1:2 with points pt 7 ps 0.5 lc 'blue' title 'Gaussian (Harmonic)', \
  'pdf_exp_harmonic.txt' using 1:2 with points pt 7 ps 0.5 lc 'red' title 'Exponential (Harmonic)', \
  f_harm(x) with lines lw 2 lc 'black' title 'Theoretical Harmonic', \
  'pdf_gauss_quartic.txt' using 1:2 with points pt 7 ps 0.5 lc 'cyan' title 'Gaussian (Quartic)', \
  'pdf_exp_quartic.txt' using 1:2 with points pt 7 ps 0.5 lc 'orange' title 'Exponential (Quartic)', \
  f_quart(x) with lines lw 2 lc 'gray' title 'Theoretical Quartic'
