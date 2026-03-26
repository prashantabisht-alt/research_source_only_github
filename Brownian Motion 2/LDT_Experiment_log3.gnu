set logscale y
set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Position Distribution: Gaussian vs Exponential Noise (t=10)'
set xrange [-30:30] # Capture tails
f(x) = exp(-x**2/(4*1.0*10.0)) / sqrt(4*pi*1.0*10.0) # t=10
plot 'pdf_gaussian.txt' u 1:2 w points title 'Gaussian Noise', \
'pdf_exponential1.txt' u 1:2 w points title 'Exponential Noise', \
f(x) w lines title 'Theoretical Gaussian'
