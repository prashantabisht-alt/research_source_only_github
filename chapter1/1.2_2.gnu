set term qt size 1100,650
set title "Velocity PDF (Maxwell–Boltzmann)"
set xlabel "v"
set ylabel "P(v)"
set grid back lw 1
set border lw 1.6
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 800

# Styles
set style line 1 lc rgb "#d62728" pt 5 ps 0.8 lw 1.4   # MC: red squares
set style line 2 lc rgb "black"   lw 2.6               # Theory

# --- Main plot ---
set multiplot
plot \
  "velocity_distributions20.txt" u 1:2 w lp ls 1 t "MC", \
  "velocity_distributions20.txt" u 1:3 w l  ls 2 t "MB theory"

# capture ranges chosen for the main panel
xL = GPVAL_X_MIN;  xR = GPVAL_X_MAX
yL = GPVAL_Y_MIN;  yR = GPVAL_Y_MAX

# --- Inset (log-y) with SAME ranges ---
set origin 0.63,0.50
set size   0.30,0.30
set logscale y
set xrange [xL:xR]
set yrange [ (yL>0 ? yL : 1e-6) : yR ]

set border lw 1
set tics out
set xlabel "v" font ",8"
set ylabel "P(v)" font ",8"
set key off

plot \
  "velocity_distributions20.txt" u 1:2 w lp ls 1 notitle, \
  "velocity_distributions20.txt" u 1:3 w l  ls 2 notitle

unset multiplot
