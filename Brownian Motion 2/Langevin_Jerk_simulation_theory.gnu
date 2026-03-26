# ===================================================================
# Gnuplot Script to Plot the Distributions from the
#           "Jerky" Langevin Simulation and Compare with Theory
#
# This script uses multiplot to generate a single output window
# containing all three distributions, each with its theoretical curve.
# ===================================================================

# --- Physical Parameters & Constants ---
PI = acos(-1.0)

# --- Begin Multiplot Layout ---
# Arrange 3 plots in 3 rows and 1 column, with a main title.
set multiplot layout 3, 1 title "Distributions for a 'Jerky' Langevin Particle vs. Theory" font ",16"

# --- Plot 1: Acceleration Distribution ---
set title "Final Acceleration Distribution P(a)"
set xlabel "Acceleration (a)"
set ylabel "Probability Density P(a)"
set key top right
set grid

# Theoretical curve for acceleration (variance = 1.0)
variance_a = 1.0
P_theory_a(x) = (1.0 / sqrt(2.0 * PI * variance_a)) * exp(-x**2 / (2.0 * variance_a))

plot 'acceleration_dist_jerky.txt' using 1:2 with points pt 7 ps 0.5 lc 'blue' title 'Simulated P(a)', \
     P_theory_a(x) with lines dashtype 2 lw 2 lc 'black' title 'Theoretical P(a)'


# --- Plot 2: Velocity Distribution ---
set title "Final Velocity Distribution P(v)"
set xlabel "Velocity (v)"
set ylabel "Probability Density P(v)"

# Theoretical curve for velocity (variance = 192 at t=50)
variance_v = 192.0
P_theory_v(x) = (1.0 / sqrt(2.0 * PI * variance_v)) * exp(-x**2 / (2.0 * variance_v))

plot 'velocity_dist_jerky.txt' using 1:2 with points pt 7 ps 0.5 lc 'red' title 'Simulated P(v)', \
     P_theory_v(x) with lines dashtype 2 lw 2 lc 'black' title 'Theoretical P(v)'


# --- Plot 3: Position Distribution ---
set title "Final Position Distribution P(x)"
set xlabel "Position (x)"
set ylabel "Probability Density P(x)"

# Theoretical curve for position (variance = 1.67e4 at t=50)
variance_x = 1.67e4
P_theory_x(x) = (1.0 / sqrt(2.0 * PI * variance_x)) * exp(-x**2 / (2.0 * variance_x))

plot 'position_dist_jerky.txt' using 1:2 with points pt 7 ps 0.5 lc 'green' title 'Simulated P(x)', \
     P_theory_x(x) with lines dashtype 2 lw 2 lc 'black' title 'Theoretical P(x)'


# --- End Multiplot ---
unset multiplot
# The plot will now appear in an interactive window.
