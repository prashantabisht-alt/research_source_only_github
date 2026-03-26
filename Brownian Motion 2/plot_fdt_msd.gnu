kBT = 1.0
m = 1.0
gamma = 0.5
D = kBT / gamma  # 2.0 for your parameters

# ---- Plot MSD ----
set xlabel 'Time t (s)'
set ylabel 'MSD ⟨x²(t)⟩ (units²)'
set title 'MSD for Full Langevin (Ballistic to Diffusive)'
set key bottom right
set grid
f_ballistic(x) = (kBT/m) * x**2
f_diffusive(x) = 2.0 * D * x

plot 'msd_gauss_fdt_full.txt' using 1:2 with points pt 7 ps 0.6 lc rgb 'blue' title 'Gaussian Noise', \
     'msd_exp_fdt_full.txt' using 1:2 with points pt 5 ps 0.6 lc rgb 'red' title 'Exponential Noise', \
     f_ballistic(x) with lines lw 2 lc rgb 'black' title 'Ballistic t²', \
     f_diffusive(x) with lines lw 2 lc rgb 'gray' title 'Diffusive 2Dt'
