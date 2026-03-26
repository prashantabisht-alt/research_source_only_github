set title "Large Deviation Function — Gaussian vs Laplace Noise"
set xlabel "Position x"
set ylabel "I(x)"
set grid
plot \
'ldf_gaussian.txt' using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'blue' title 'Gaussian', \
'ldf_exponential.txt' using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'red' title 'Laplace'
