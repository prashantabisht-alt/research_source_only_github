# ===================== Parameters (match your sim) =====================
D    = 1.0
dt   = 0.001
Tmax = 200
T    = Tmax*dt                 # = 0.2
pi   = 3.141592653589793
a    = sqrt(2.0*D*dt)          # a = sqrt(2 D dt)

# ===================== Helpers for empirical I_T =======================
Y(colx)  = (colx)/T
IT(colp) = (colp>0 ? (-1.0/T)*log(colp) : 1/0)

# ===================== Gaussian theory (finite-T) ======================
# Exact finite-T normalization constant:
Cgauss(T,D) = (1.0/(2.0*T))*log(4.0*pi*D*T)
Igauss(y)   = y*y/(4.0*D) + Cgauss(T,D)

# ===================== Exponential-increment theory ====================
# Asymptotic LDP rate (T→∞), mapped to y = x/T; domain y > -a/dt
Iexp_asym(y) = (y > -a/dt ? (1.0/dt) * ( (dt/a)*y - log(1.0 + (dt/a)*y) ) : 1/0)

# Add a constant to align the center (like Gaussian). You can keep 2.3 or auto-fit below.
Cexp = 2.3
Iexp_const(y) = (y > -a/dt ? Iexp_asym(y) + Cexp : 1/0)

# --- Finite-T saddlepoint correction (improves tail agreement):
# SCGF: lambda(k) = [-a k - ln(1 - a k)] / dt
# Stationary point: k*(y) = y / (2D + a y)
# lambda''(k*) = (2D + a y)^2 / (2D)
Iexp_corr(y) = (y > -a/dt ? Iexp_asym(y) + (1.0/(2.0*T))*log( 2.0*pi*T * ((2.0*D + a*y)**2 / (2.0*D)) ) : 1/0)

# ===================== (Optional) auto-fit Cexp near center =============
# Un-comment to let gnuplot find best constant (recommended).
# y0 = 2.0
# fit (y > -a/dt ? Iexp_asym(y) + Cexp : 1/0) \
#     "exponential_shifted.txt" \
#     u (Y($1)):( ($2>0 && abs(Y($1))<=y0) ? IT($2) : 1/0 ) via Cexp
# print sprintf("Fitted Cexp = %.6f over |y|<=%.2f", Cexp, y0)

# ===================== Plot setup ======================================
set term qt size 1100,650
set title sprintf("Scaled rate function I_T(y) at T = %.3g", T)
set xlabel "y = x / T"
set ylabel "I_T(y) = -(1/T) ln P_T(x = y T)"
set grid back lw 1
set border lw 1.8
set tics out
set key left box opaque width -2 spacing 1.2
set samples 800
set xrange [-13:13]
set yrange [0:45]

# Styles
set style line 1 lc rgb "#1f77b4" lw 2.2              # Gaussian (emp)
set style line 2 lc rgb "#d62728" lw 2.2              # Exp-shifted (emp)
set style line 4 lc rgb "#000000" lw 2.8 dt 2         # Gaussian theory (dashed)
set style line 5 lc rgb "#2ca02c" lw 2.8              # Exp theory (asym + const)
set style line 6 lc rgb "#228833" lw 2.8 dt 3         # Exp theory (finite-T corrected)

# ===================== Plot ============================================
plot \
  "gaussian_simulation1.txt"  u (Y($1)):(IT($2)) w l ls 1 t "Gaussian (emp)", \
  "exponential_shifted.txt"   u (Y($1)):(IT($2)) w l ls 2 t "Exp-shifted (emp)", \
  Igauss(x)                   w l ls 4 t "Gaussian theory", \
  Iexp_const(x)               w l ls 5 t "Exp theory (asym + const)", \
  Iexp_corr(x)                w l ls 6 t "Exp theory (finite-T corr)"


