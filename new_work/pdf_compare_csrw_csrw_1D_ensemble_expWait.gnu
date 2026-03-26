# ==== Files ====
base_tag_ctrw = "ctrw_1D_ensemble_expWait_exactT_lambda1p0_D1p0_tmax10"
base_tag_csrw = "csrw_1D_ensemble_exactT_lambdaNA_D1p0_tmax10"

f_ctrw = "xPDF_".base_tag_ctrw.".txt"
f_csrw = "xPDF_".base_tag_csrw.".txt"

# ==== Parameters (match code) ====
D    = 1.0
tmax = 10000.0
pi   = 3.141592653589793

# Theory: N(0, 2 D tmax)
sigma2 = 2.0*D*tmax
sigma  = sqrt(sigma2)
g(x)   = 1.0/(sqrt(2*pi)*sigma) * exp(-0.5*(x*x)/sigma2)

# ==== Styling ====
set term qt size 1100,650
set title sprintf("Snapshot P(x, t=%.0f): CTRW vs CSRW vs N(0,2Dt)", tmax)
set xlabel "x"
set ylabel "P(x, t)"
set grid back lw 1
set border lw 1.5
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 600

set style line 1 lc rgb "#d62728" pt 5 ps 0.9 lw 1.4   # CTRW: red squares + thin line
set style line 2 lc rgb "#1f77b4" pt 7 ps 0.9 lw 1.4   # CSRW: blue triangles + thin line
set style line 3 lc rgb "black"   lw 2.8               # Theory: solid black

plot \
  f_ctrw u 1:2 w lp ls 1 t "CTRW (exp waits, exact T)", \
  f_csrw u 1:2 w lp ls 2 t "CSRW (discrete-time exact T)", \
  g(x)   w  l  ls 3 t "Gaussian theory"
