
set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Position Distribution: Gaussian vs Exponential Noise'

# Theoretical Gaussian for D=1, t=1
f(x) = exp(-x**2/(4.0*1.0*1.0)) / sqrt(4.0*pi*1.0*1.0)

plot 'pdf_gaussian.txt' using 1:2 with points pt 2 ps 0.8 lc rgb 'blue' title 'Gaussian Noise', \
     'pdf_exponential.txt' using 1:2 with points pt 2 ps 0.8 lc rgb 'red' title 'Exponential Noise', \
     f(x) with lines lw 2 lc rgb 'black' title 'Theoretical Gaussian'
