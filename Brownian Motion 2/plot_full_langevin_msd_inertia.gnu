set xlabel 'Time t (s)'
set ylabel 'MSD <x^2> (units^2)'
set title 'MSD for Full Langevin with Inertia Only (t = 50)'
set key bottom right

# Theoretical curves
kBT = 1.0
m = 1.0
f_ballistic(x) = (kBT / m) * x**2      # Ballistic regime
f_diffusive(x) = 2.0 * x               # Diffusive regime, D = kBT / gamma = 2.0

# Plot: points for data, lines for theory
plot 'msd_gauss_inertia_only.txt' using 1:2 with points pt 7 ps 0.5 lc 'blue' title 'Gaussian Noise', \
     'msd_exp_inertia_only.txt' using 1:2 with points pt 5 ps 0.5 lc 'red' title 'Exponential Noise', \
     f_ballistic(x) with lines lw 2 lc 'black' title 'Ballistic t^2', \
     f_diffusive(x) with lines lw 2 lc 'gray' title 'Diffusive 2Dt'
