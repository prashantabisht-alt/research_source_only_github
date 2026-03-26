set title "OU-Colored Noise (Overdamped) — Gaussian vs Laplace, Different τ_c"
set xlabel 'Position x'
set ylabel 'P(x)'
set grid
set key top right

# linepoints style: pt 2 = cross, ps 0.6, lw 1 = thin line
plot \
'pdf_gauss_tau0.txt'    using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'blue'  title 'Gaussian τ_c=0', \
'pdf_gauss_tau0.1.txt'  using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'cyan'  title 'Gaussian τ_c=0.1', \
'pdf_gauss_tau1.0.txt'  using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'navy'  title 'Gaussian τ_c=1.0', \
'pdf_exp_tau0.txt'      using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'red'   title 'Laplace τ_c=0', \
'pdf_exp_tau0.1.txt'    using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'orange' title 'Laplace τ_c=0.1', \
'pdf_exp_tau1.0.txt'    using 1:2 with linespoints pt 2 ps 0.6 lw 1 lc rgb 'dark-red' title 'Laplace τ_c=1.0'
