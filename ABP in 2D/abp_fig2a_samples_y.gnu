reset
set term qt size 1100,650

set xlabel "y / σ_y"
set ylabel "P(y,t) σ_y"
set key left box opaque
set grid back lw 1

# Standard normal for reference
normal(x) = 1.0/sqrt(2.0*pi)*exp(-0.5*x**2)

set multiplot

# --- Main plot (linear y) ---
plot "y_scaled_t0.25.txt" u 1:2 w lp t "t=0.25", \
     "y_scaled_t0.50.txt" u 1:2 w lp t "t=0.50", \
     "y_scaled_t1.00.txt" u 1:2 w lp t "t=1.00", \
     normal(x)             w l  lw 2 t "N(0,1)"

# capture the auto ranges chosen for the main panel
xL = GPVAL_X_MIN;  xR = GPVAL_X_MAX
yL = GPVAL_Y_MIN;  yR = GPVAL_Y_MAX

# --- Inset (log-y) with same numeric ranges ---
set origin 0.60,0.52
set size   0.36,0.36
set logscale y
set xrange [xL:xR]
# ensure positive lower bound for log-scale
set yrange [ (yL>0 ? yL : 1e-8) : yR ]
set border lw 1
set tics out
set xlabel "y/σ_y" font ",8"
set ylabel "Pσ_y"  font ",8"
set key off

plot "y_scaled_t0.25.txt" u 1:2 w lp, \
     "y_scaled_t0.50.txt" u 1:2 w lp, \
     "y_scaled_t1.00.txt" u 1:2 w lp, \
     normal(x)             w l  lw 2

unset multiplot
