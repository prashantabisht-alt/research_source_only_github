set xlabel 'Time t (s)'
set ylabel 'MSD ⟨x²(t)⟩ (units²)'
set title 'MSD for Overdamped Langevin'
set key bottom right
f_msd(x) = 2.0 * x  # Theoretical: 2 D t with D = kBT / γ = 1.0
plot 'msd_gauss_fdt_overdamped.txt' using 1:2 with points pointtype 7 pointsize 1 lc rgb 'blue' title 'Gaussian Noise', \
     'msd_exp_fdt_overdamped.txt' using 1:2 with points pointtype 5 pointsize 1 lc rgb 'red' title 'Exponential Noise', \
     f_msd(x) with lines lw 2 lc rgb 'black' title 'Theoretical 2Dt'

# --- Second subplot: Autocorrelation ---
set xlabel 'Time t (s)'
set ylabel 'Autocorrelation ⟨x(t)x(0)⟩ (units²)'
set title 'Position Autocorrelation (Overdamped Limit)'
f_autocorr(x) = 1.0  # Constant correlation at t=0 (for normalized units)
plot 'autocorr_gauss_fdt_overdamped.txt' using 1:2 with points pointtype 7 pointsize 1 lc rgb 'blue' title 'Gaussian Noise', \
     'autocorr_exp_fdt_overdamped.txt' using 1:2 with points pointtype 5 pointsize 1 lc rgb 'red' title 'Exponential Noise', \
     f_autocorr(x) with lines lw 2 lc rgb 'black' title 'Theoretical (Constant)'

unset multiplot
