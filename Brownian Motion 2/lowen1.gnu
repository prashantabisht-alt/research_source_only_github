set multiplot layout 1,2 title "AOUP u(t) statistics"

# --- Autocorrelation C(tau) ---
set title "Velocity autocorrelation C(τ)"
set xlabel "τ"
set ylabel "C(τ)"
set grid
set key top right
plot \
  "u_autocorr.txt" using 1:2 with lines lc rgb "blue" lw 2 title "C_est", \
  "u_autocorr.txt" using 1:3 with lines lc rgb "red"  lw 2 dt 2 title "C_theory"

# --- Memory kernel M(tau) ---
set title "Memory kernel M(τ) = γ² C(τ)"
set xlabel "τ"
set ylabel "M(τ)"
set grid
set key top right
plot \
  "M_kernel.txt" using 1:2 with lines lc rgb "green" lw 2 title "M_est", \
  "M_kernel.txt" using 1:3 with lines lc rgb "orange" lw 2 dt 2 title "M_theory"

unset multiplot
