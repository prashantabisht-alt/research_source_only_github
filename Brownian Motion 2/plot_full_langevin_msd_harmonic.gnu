set xlabel 'Time t (s)'
set ylabel 'MSD <x^2> (units^2)'
set title 'MSD in Harmonic Potential (Full Langevin, x = 50)'
set key top left
set grid

# Optional theoretical curve for harmonic oscillator MSD:
f_harm(x) = 5.0 * (1.0 - exp(-x / 2.0))  # tau = m/gamma = 2.0, k = 0.2, so saturation ≈ 5

# Plotting only the harmonic data
plot 'msd_harm_gauss.txt' using 1:2 with points pt 7 ps 0.7 lc rgb 'blue' title 'Harmonic (Gaussian)', \
     'msd_harm_exp.txt'   using 1:2 with points pt 5 ps 0.7 lc rgb 'red' title 'Harmonic (Exponential)', \
     f_harm(x) with lines lw 2 lc rgb 'black' title 'Theoretical Harmonic'
