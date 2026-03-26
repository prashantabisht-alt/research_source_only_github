# --- config: adjust file name to match your base_tag ---
pdf  = "xPDF_ctrw_1D_ensemble_expWait_exactT_lambda1p0_D1p0_tmax10.txt"

# parameters (must match your simulation)
D     = 1.0
tmax  = 10000.0      # you set in code
pi    = 3.141592653589793

# theory: Gaussian N(0, 2 D tmax)
sigma2 = 2.0*D*tmax
sigma  = sqrt(sigma2)
g(x)   = 1.0/(sqrt(2*pi)*sigma) * exp(-0.5*(x/sigma)**2)

# --- plot PDF vs Gaussian ---
set term qt 0
set title sprintf("Snapshot P(x, tmax=%.0f) vs N(0, 2Dt)", tmax)
set xlabel "x"
set ylabel "P(x, tmax)"
set grid
plot pdf u 1:2 w lp lw 2 t "Monte Carlo PDF", \
     g(x) w l dt 2 lw 2 t sprintf("Gaussian theory, σ^2=%.0f", sigma2)
