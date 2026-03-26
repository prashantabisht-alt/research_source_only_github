### ---- USER PARAMS: match to your run ----
file  = "xPDF_ctrw_1D_ensemble_expWait_exactT_lambda1p0_D1p0_tmax10.txt"  # update name if needed
D     = 1.0
tmax  = 10000.0
nbins = 200
xmin  = -800.0
xmax  =  800.0

### ---- Derived ----
dx    = (xmax - xmin)/nbins
g(x)  = 1.0/sqrt(4.0*pi*D*tmax) * exp(-x*x/(4.0*D*tmax))

### ---- Stats from simulation ----
# file: col1 = x center, col2 = PDF value
stats file using ($1*$2)    name "M1" nooutput
stats file using ($1*$1*$2) name "M2" nooutput
mean_x = M1_sum * dx
var_x  = M2_sum * dx - mean_x**2

print sprintf("Empirical mean  = %.6f   (theory: 0)", mean_x)
print sprintf("Empirical var   = %.6f   (theory: %.6f = 2*D*tmax)", var_x, 2*D*tmax)

### ---- Plot ----
set term qt 0
set title sprintf("P(x,t) at t=%.0f  (D=%.1f, nbins=%d)", tmax, D, nbins)
set xlabel "x"
set ylabel "P(x,t)"
set grid
set xrange [xmin:xmax]

plot file using 1:2 with points lw 2 title "Simulation PDF", \
     g(x) with lines lw 2 title "Theory Gaussian"
