set logscale y
set xlabel 'Position x'
set ylabel 'P(x)'
set title 'PDF — Gaussian vs Laplace Noise (Overdamped)'
f_theory(x) = exp(-x**2/(4*1.0*100.0)) / sqrt(4*pi*1.0*100.0)   # D=1, t=100
plot 'pdf_gaussian.txt' using 1:2 with points pt 7 ps 0.5 lc 'blue' title 'Gaussian', \
     'pdf_laplace.txt' using 1:2 with points pt 7 ps 0.5 lc 'red' title 'Laplace', \
     f_theory(x) with lines lw 2 lc 'black' title 'Theoretical Gaussian'

set output 'ldf_plot.png'
unset logscale y
set xlabel 'Position x'
set ylabel 'I(x)'
set title 'Large Deviation Function — Gaussian vs Laplace'
plot 'ldf_gaussian.txt' using 1:2 with linespoints pt 7 ps 0.5 lc 'blue' title 'Gaussian', \
     'ldf_laplace.txt' using 1:2 with linespoints pt 7 ps 0.5 lc 'red' title 'Laplace'
