
kBT = 1.0
m = 1.0
gamma = 0.5
tau = m / gamma  # 2.0 for your parameters

# ---- Plot VACF ----
set xlabel 'Time t (s)'
set ylabel 'VACF ⟨v(t)v(0)⟩ (units²)'
set title 'Velocity Autocorrelation for Full Langevin'
set key top right
set grid
f_autocorr(x) = (kBT/m) * exp(-x / tau)

plot 'autocorr_gauss_fdt_full.txt' using 1:2 with points pt 7 ps 0.6 lc rgb 'blue' title 'Gaussian Noise', \
     'autocorr_exp_fdt_full.txt' using 1:2 with points pt 5 ps 0.6 lc rgb 'red' title 'Exponential Noise', \
     f_autocorr(x) with lines lw 2 lc rgb 'black' title 'Theoretical exp(-t/τ)'
