# ===== Brownian motion: PDF at tmax with log-y inset =====
fdata = "gaussian_simulation1.txt"

# --- simulation params (match your Fortran) ---
dt    = 0.001
Tmax  = 200
D     = 1.0
tmax  = dt*Tmax
pi    = 3.141592653589793

# --- Gaussian theory: N(0, 2 D tmax) ---
sigma2 = 2.0*D*tmax
sigma  = sqrt(sigma2)
g(x)   = 1.0/(sqrt(2*pi)*sigma) * exp(-0.5*(x*x)/sigma2)

# --- aesthetics ---
set term qt size 1100,650
set title sprintf("Brownian PDF at t = %.3f  (D=%.1f):  MC vs N(0, 2Dt)", tmax, D)
set xlabel "x"
set ylabel "P(x, t)"
set grid back lw 1
set border lw 1.6
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 800

# styles
set style line 1 lc rgb "red" pt 7 ps 0.6 lw 1.4   # MC: blue triangles
set style line 2 lc rgb "black"   lw 2.8               # Theory: solid black

# --- Main plot ---
set multiplot

plot \
  fdata u 1:2 w lp ls 1 t "MC (histogram → PDF)", \
  g(x)        w  l ls 2 t sprintf("Gaussian theory (σ^2 = %.3f)", sigma2)

# --- Inset (log-scale y) ---
set origin 0.55,0.45     # position of inset (0–1 screen coords)
set size   0.4,0.4       # size of inset
set logscale y
unset key
unset title
set xlabel ""
set ylabel ""
set grid back lw 0.8
set tics out

replot

unset multiplot
