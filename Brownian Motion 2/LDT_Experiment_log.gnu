# Gnuplot Script: LDT_Experiment_log.gnu
# This script plots data on a log scale in an interactive window.

# --- Plot Aesthetics ---
set title "Probability Distribution (Log Scale)"
set xlabel "Final Position (x)"
set ylabel "Probability Density P(x)"
set grid
set key top left

# --- KEY COMMAND: Set Y-axis to Logarithmic Scale ---
set logscale y
set format y "10^{%L}" # Optional: makes y-axis labels look nice

# --- Define Constants and Theoretical Functions ---
N_steps = 50.0
diff = 1.0
dt = 0.001
sigma_sq = N_steps * 2.0 * diff * dt
gauss(x) = (1.0 / sqrt(2*pi*sigma_sq)) * exp(-x**2 / (2.0*sigma_sq))
lambda = 1.0 / sqrt(2.0 * diff * dt)
k = N_steps
erlang(x) = (x > 0) ? (lambda**k / gamma(k)) * x**(k-1) * exp(-lambda*x) : 1/0

# --- Plot Command ---
# This plots the data from your .txt files directly to the interactive window
plot 'pdf_gaussian.txt' using 1:2 with lines title "Simulation (Gaussian)", \
     gauss(x) with lines lw 2 lc 'black' title "Theoretical Gaussian", \
     'pdf_exponential.txt' using 1:2 with lines title "Simulation (Exponential)", \
     erlang(x) with lines lw 2 dt '--' lc 'black' title "Theoretical Erlang"

print "Plot displayed. Close the plot window or type 'exit' to quit."
