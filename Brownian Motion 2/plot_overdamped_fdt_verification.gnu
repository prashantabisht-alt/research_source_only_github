
# --- Plot Aesthetics ---
set title "FDT Verification for Overdamped Motion: <x^2> = 2Dt" font ",16"
set xlabel "Time (t)"
set ylabel "Mean Squared Displacement <x^2>"
set key top left
set grid

# --- Define Theoretical Parameters (must match Fortran code) ---
kBT = 1.0
gamma = 1.0
D = kBT / gamma
theoretical_line(x) = 2.0 * D * x

# --- Plot Command ---
plot 'msd_overdamped_gaussian.txt' using 1:2 with points pt 7 ps 0.5 lc 'blue' title 'Simulation (Gaussian Noise)', \
     'msd_overdamped_exponential.txt' using 1:2 with points pt 5 ps 0.5 lc 'red' title 'Simulation (Exponential Noise)', \
     theoretical_line(x) with lines dashtype 2 lw 2 lc 'black' title sprintf("Theoretical Line (2Dt, D=%.1f)", D)
