# --- config: adjust file name to match your base_tag ---
pdf  = "xPDF_ctrw_1D_ensemble_expWait_exactT_lambda1p0_D1p0_tmax10.txt"

# parameters (must match your simulation)
D     = 1.0
tmax  = 10000.0
pi    = 3.141592653589793

# theory: Gaussian N(0, 2 D tmax)
sigma2 = 2.0*D*tmax
sigma  = sqrt(sigma2)
g(x)   = 1.0/(sqrt(2*pi)*sigma) * exp(-0.5*(x/sigma)**2)

# --- aesthetics ---
set term qt size 1100,650
set title sprintf("Snapshot P(x, t=%.0f) vs N(0, 2Dt)", tmax)
set xlabel "x"
set ylabel "P(x, t)"
set grid back lw 1
set border lw 1.5
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 600   # smooth function curve

# styles (color-blind friendly)
set style line 1 lc rgb "#d62728" pt 5 ps 0.9 lw 1.4  # MC: red squares + thin line
set style line 2 lc rgb "black"   lw 2.8              # Theory: solid black

set multiplot

# --- main plot (xrange fixed) ---
set xrange [-600:800]
plot \
  pdf u 1:2 w lp ls 1 t "Monte Carlo PDF", \
  g(x)    w  l ls 2 t sprintf("Gaussian theory (σ^2 = %.0f)", sigma2)

# --- inset (log-y, auto x) ---
set autoscale x   # <-- reset xrange to auto
set origin 0.58,0.50
set size   0.36,0.36
set logscale y
set border lw 1
set tics out
set xlabel "x" font ",8"
set ylabel "P" font ",8"
set key off

plot \
  pdf u 1:2 w lp ls 1, \
  g(x)    w  l ls 2

unset multiplot
