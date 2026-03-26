set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Position Distribution: Gaussian vs laplace Noise'
# Theoretical Gaussian: P(x) = 1/sqrt(4 pi D t) * exp(-x^2/(4 D t)), D=1, t=Tmax*dt=200*0.001=0.2
f(x) = exp(-x**2/(4*1.0*0.2)) / sqrt(4*pi*1.0*0.2)
plot 'gaussian_simulation1.txt' using 1:2 with points title 'Gaussian Noise', \
     'laplace.txt' using 1:2 with points title 'laplace Noise', \
     f(x) with lines title 'Theoretical Gaussian'
