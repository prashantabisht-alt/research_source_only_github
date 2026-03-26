set logscale y
set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Overdamped Position Distribution with Harmonic Potential (t=10)'
set xrange [-10:10]
set key bottom right
f(x) = exp(-0.5 * x**2) / sqrt(2.0 * 3.1415926535)  # Theoretical P(x) ∝ exp(-k x^2 / 2 k_B T), k=1
plot 'pdf_gauss_harmonic.txt' using 1:2 with points title 'Gaussian Noise', \
     'pdf_exp_harmonic.txt' using 1:2 with points title 'Exponential Noise', \
     f(x) with lines title 'Theoretical Boltzmann'
