# === Brownian final-time PDFs: data + Gaussian theory ===
# Match these to your simulation params:
D    = 1.0
dt   = 0.001
Tmax = 200
tfin = dt*Tmax

pi = 3.141592653589793
sigma2 = 2.0*D*tfin
sigma  = sqrt(sigma2)
g(x)   = 1.0/(sqrt(2*pi)*sigma) * exp(-0.5*(x*x)/sigma2)

set term qt size 1100,650
set title sprintf("Final-time P(x, T=%.3g) — Gaussian vs Exponential (theory σ^2=%.3g)", tfin, sigma2)
set xlabel "x"
set ylabel "P(x, T)"
set grid back lw 1
set border lw 1.8
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 600

# Styles
set style line 1 lc rgb "#1f77b4" pt 7 ps 0.8 lw 2    # Gaussian increments (MC)
set style line 2 lc rgb "#d62728" pt 5 ps 0.8 lw 2    # Exponential (shifted+scaled) (MC)
set style line 4 lc rgb "black"   lw 2.8              # Theory: solid black

plot \
  "gaussian_simulation1.txt" u 1:2 w lp ls 1 t "Gaussian (MC)", \
  "exponential_shifted.txt"  u 1:2 w lp ls 2 t "Exponential (MC)", \
  g(x) w l ls 4 t "Gaussian theory  N(0, 2Dt)", \
  "exp_theory.txt" u 1:2 w l lw 2 lc rgb "#228833" t "Exponential theory (exact finite-T)"

