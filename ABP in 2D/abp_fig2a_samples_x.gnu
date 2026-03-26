reset
set term qt size 1100,650

# --- Main plot ---
set xlabel "a1 = (v0 t - x)/(v0 DR t^2)"
set ylabel "P(x,t) v0 DR t^2"
set key left box opaque
set grid back lw 1

set multiplot

# main panel with xrange capped at 1.5
set xrange [:1.5]
plot "a1_edge_t0.25.txt" u 1:2 w lp t "t=0.25", \
     "a1_edge_t0.50.txt" u 1:2 w lp t "t=0.50", \
     "a1_edge_t1.00.txt" u 1:2 w lp t "t=1.00", \
     "fx_theory.txt"     u 1:2 w l  lw 2 t "theory f_x"

# --- Inset (log-y) with SAME x-range ---
xL = GPVAL_X_MIN;  xR = GPVAL_X_MAX
yL = GPVAL_Y_MIN;  yR = GPVAL_Y_MAX

set origin 0.58,0.50   # top-right corner
set size   0.38,0.38
set logscale y
set xrange [xL:xR]
set yrange [ (yL>0 ? yL : 1e-6) : yR ]
set border lw 1
set tics out
set xlabel "a1" font ",8"
set ylabel "P"  font ",8"
set key off

plot "a1_edge_t0.25.txt" u 1:2 w lp ls 1, \
     "a1_edge_t0.50.txt" u 1:2 w lp ls 2, \
     "a1_edge_t1.00.txt" u 1:2 w lp ls 3, \
     "fx_theory.txt"     u 1:2 w l  lw 2

unset multiplot
