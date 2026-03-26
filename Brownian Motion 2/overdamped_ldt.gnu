set title "Overdamped LDT — Gaussian vs Laplace Noise"
set xlabel "Position x"
set ylabel "P(x)"
set logscale y
set grid
set key top right

# Parameters
D = 1.0
tmax = 1000 * 0.01   # Tmax * dt from code
f_theory(x) = exp(-x**2/(4*D*tmax)) / sqrt(4*pi*D*tmax)

plot \
'pdf_gaussian.txt' using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'blue' title 'Gaussian', \
'pdf_exponential.txt' using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'red'  title 'Laplace', \
f_theory(x) with lines lw 2 lc rgb 'black' title 'Theoretical Gaussian'
