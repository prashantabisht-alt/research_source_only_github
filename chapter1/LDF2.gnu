# ===== Parameters (match your sim) =====
D    = 1.0
dt   = 0.001
Tmax = 200
T    = Tmax*dt       # = 0.2


Igauss(y) = y*y/(4.0*D) + 2.3
Y(colx)   = (colx)/T
IT(colp)  = (colp>0 ? (-1.0/T)*log(colp) : 1/0)

set term qt size 1100,650
set title sprintf("Scaled rate function I_T(y) at T = %.3g", T)
set xlabel "y = x / T"
set ylabel "I_T(y) = -(1/T) ln P_T(x = y T)"
set grid back lw 1
set border lw 1.8
set tics out
set key left box opaque width -2 spacing 1.2
set samples 600

set style line 1 lc rgb "#1f77b4"  lw 2.2                # Gaussian (emp)
set style line 2 lc rgb "#d62728"  lw 2.2                # Exp-shifted (emp)
set style line 4 lc rgb "black"    lw 2.8 dt 2           # Gaussian theory (dashed)

plot \
  "gaussian_simulation1.txt"  u (Y($1)):(IT($2)) w l ls 1 t "Gaussian (emp)", \
  "exponential_shifted.txt"   u (Y($1)):(IT($2)) w l ls 2 t "Exp-shifted (emp)", \
  Igauss(x)                   w l ls 4 t "Gaussian theory"
