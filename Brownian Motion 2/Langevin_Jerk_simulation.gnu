# ===================================================================
# Gnuplot Script to Plot the Distributions from the
#           "Jerky" Langevin Simulation
#
# This script uses multiplot to generate a single output window
# containing all three distributions arranged vertically.
# ===================================================================

# --- Physical Parameters (must match Fortran code) ---
m = 1.0
kBT = 1.0
PI = acos(-1.0)

# --- Begin Multiplot Layout ---
# Arrange 3 plots in 3 rows and 1 column, with a main title.
set multiplot layout 3, 1 title "Distributions for a 'Jerky' Langevin Particle (m*a_dot + gamma*a = xi)" font ",16"

# --- Plot 1: Acceleration Distribution ---
set title "Final Acceleration Distribution P(a)"
set xlabel "Acceleration (a)"
set ylabel "Probability Density P(a)"
set key top right
set grid
# set logscale y  -- REMOVED as requested

# Theoretical curve for acceleration (Ornstein-Uhlenbeck process)
variance_a_theory = kBT / m
norm_factor_a = 1.0 / sqrt(2.0 * PI * variance_a_theory)
P_theory_a(x) = norm_factor_a * exp(-x**2 / (2.0 * variance_a_theory))

plot 'acceleration_dist_jerky.txt' using 1:2 with points pt 7 ps 0.5 lc 'blue' title 'Simulated P(a)', \
     P_theory_a(x) with lines dashtype 2 lw 2 lc 'black' title 'Theoretical P(a) (Gaussian)'


# --- Plot 2: Velocity Distribution ---
set title "Final Velocity Distribution P(v)"
set xlabel "Velocity (v)"
set ylabel "Probability Density P(v)"
# Other aesthetics are inherited from the previous plot

plot 'velocity_dist_jerky.txt' using 1:2 with points pt 7 ps 0.5 lc 'red' title 'Simulated P(v)'


# --- Plot 3: Position Distribution ---
set title "Final Position Distribution P(x)"
set xlabel "Position (x)"
set ylabel "Probability Density P(x)"

plot 'position_dist_jerky.txt' using 1:2 with points pt 7 ps 0.5 lc 'green' title 'Simulated P(x)'


# --- End Multiplot ---
unset multiplot
# The plot will now appear in an interactive window.

