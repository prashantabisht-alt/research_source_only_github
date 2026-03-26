set term qt size 1100,650
set title "Position PDF at t=50 (D=2)"

set xlabel "x"
set ylabel "P(x,t)"
set grid back lw 1
set border lw 1.6
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 800

# Styles
set style line 3 lc rgb "#1f77b4" pt 7 ps 0.8 lw 1.4   # MC
set style line 4 lc rgb "black"   lw 2.6               # Theory

# --- Main plot ---
set multiplot
plot \
  "position_distributions20.txt" u 1:2 w lp ls 3 t "MC", \
  "position_distributions20.txt" u 1:3 w l  ls 4 t "Gaussian"

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
set xlabel "x" font ",8"
set ylabel "P(x)" font ",8"
set key off

plot \
  "position_distributions20.txt" u 1:2 w lp ls 3 notitle, \
  "position_distributions20.txt" u 1:3 w l  ls 4 notitle

unset multiplot
