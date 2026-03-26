set xlabel 'Time t (s)'
set ylabel 'MSD <x^2> (units^2)'
set title 'MSD for Full Langevin (FDT Verification)'
set key bottom right
f_msd(x) = 2.0 * x  # Long-time 2 D t, D = kBT / gamma = 2.0
plot 'msd_gauss_fdt_full.txt' using 1:2 with points pt 7 ps 0.6 lc rgb 'blue' title 'Gaussian Noise', \
     'msd_exp_fdt_full.txt' using 1:2 with points pt 5 ps 0.6 lc rgb 'red' title 'Exponential Noise', \
     f_msd(x) with lines lw 2 lc rgb 'black' title 'Theoretical 2Dt'

# --- Second subplot: Velocity Autocorrelation ---
set xlabel 'Time t (s)'
set ylabel 'Autocorrelation <v(t)v(0)> (units^2)'
set title 'Velocity Autocorrelation (FDT Verification)'
f_autocorr(x) = (kBT / m) * exp(-x / tau)
plot 'autocorr_gauss_fdt_full.txt' using 1:2 with points pt 7 ps 0.6 lc rgb 'blue' title 'Gaussian Noise', \
     'autocorr_exp_fdt_full.txt' using 1:2 with points pt 5 ps 0.6 lc rgb 'red' title 'Exponential Noise', \
     f_autocorr(x) with lines lw 2 lc rgb 'black' title 'Theoretical exp(-t/tau)'

unset multiplot
