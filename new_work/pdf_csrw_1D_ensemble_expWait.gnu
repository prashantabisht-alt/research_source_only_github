# ===== CSRW-only PDF plot =====
# Adjust the base tag if you changed it in csrw_1D_ensemble_exactT.f90
base_tag_csrw = "csrw_1D_ensemble_exactT_lambdaNA_D1p0_tmax10"
f_csrw = "xPDF_".base_tag_csrw.".txt"

# Simulation params (match your code)
D    = 1.0
tmax = 10000.0
pi   = 3.141592653589793

# Gaussian theory: N(0, 2 D tmax)
sigma2 = 2.0*D*tmax
sigma  = sqrt(sigma2)
g(x)   = 1.0/(sqrt(2*pi)*sigma) * exp(-0.5*(x*x)/sigma2)

# --- aesthetics ---
set term qt size 1100,650
set title sprintf("CSRW snapshot P(x, t=%.0f) vs Gaussian N(0,2Dt)", tmax)
set xlabel "x"
set ylabel "P(x, t)"
set grid back lw 1
set border lw 1.5
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 600

# styles
set style line 1 lc rgb "#1f77b4" pt 7 ps 0.9 lw 1.4   # CSRW: blue triangles
set style line 2 lc rgb "black"   lw 2.8               # Theory: solid black

set multiplot

# --- Main plot (xrange fixed) ---
set xrange [-600:800]
plot \
  f_csrw u 1:2 w lp ls 1 t "CSRW PDF (MC)", \
  g(x)   w  l  ls 2 t "Gaussian theory"

# --- Inset (log-y, auto-x) ---
set autoscale x     # reset to auto for inset
set origin 0.55,0.45
set size   0.4,0.4
set logscale y
unset key
unset title
set xlabel "" font ",8"
set ylabel "" font ",8"
set grid back lw 0.8
set tics out

plot \
  f_csrw u 1:2 w lp ls 1, \
  g(x)   w  l  ls 2

unset multiplot
