reset
set term qt size 1100,650

set xlabel "z = (v0^2 t^2 - r^2)/(2 DR v0^2 t^3)"
set ylabel "P(r^2,t) · 2 DR v0^2 t^3"
set key left box opaque
set grid back lw 1

set multiplot

# --- Main (linear scale) ---
plot "fr_t0.25.txt" u 1:2 w lp t "t=0.25", \
     "fr_t0.50.txt" u 1:2 w lp t "t=0.50", \
     "fr_t1.00.txt" u 1:2 w lp t "t=1.00", \
     "fr_theory.txt" u 1:2 w l  lw 2 lc rgb "red" t "theory f_r(z)"

# capture auto ranges
xL = GPVAL_X_MIN;  xR = GPVAL_X_MAX
yL = GPVAL_Y_MIN;  yR = GPVAL_Y_MAX

# --- Inset (log-y) ---
set origin 0.60,0.52
set size   0.36,0.36
set logscale y
set xrange [xL:xR]
set yrange [ (yL>0 ? yL : 1e-8) : yR ]   # clamp if yL ≤ 0
set border lw 1
set tics out
set xlabel "z" font ",8"
set ylabel "P" font ",8"
set key off

plot "fr_t0.25.txt" u 1:2 w lp, \
     "fr_t0.50.txt" u 1:2 w lp, \
     "fr_t1.00.txt" u 1:2 w lp, \
     "fr_theory.txt" u 1:2 w l lw 2 lc rgb "red"

unset multiplot
