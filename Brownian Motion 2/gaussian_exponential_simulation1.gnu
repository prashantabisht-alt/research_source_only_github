
set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Position Distribution: Gaussian vs Exponential Noise'
# Theoretical Gaussian: P(x) = 1/sqrt(4 pi D t) * exp(-x^2/(4 D t)), D=1, t=Tmax*dt=1000*0.001=1
f(x) = exp(-x**2/(4*1.0*1.0)) / sqrt(4*pi*1.0*1.0)
plot 'pdf_gaussian.txt' using 1:2 with points title 'Gaussian Noise', \
     'pdf_exponential.txt' using 1:2 with points title 'Exponential Noise', \
     f(x) with lines title 'Theoretical Gaussian'
