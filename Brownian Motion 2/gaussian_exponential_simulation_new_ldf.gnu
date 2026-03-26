reset
# ---- set your run params here ----
D = 1.0                 # your 'diff'
dt = 0.001              # your time step
Tmax = 200              # your steps
T = Tmax * dt           # total physical time

# ---- axes + labels ----
set term qt
set xlabel "y = x / T"
set ylabel "I_T(y) = -(1/T) log P_T(x=y T)"
set key left
set grid

# Theoretical LDF for Gaussian: I(y) = y^2 / (4D)
Igauss(y) = y*y/(4.0*D)

# Build columns on the fly:
# 1st col in files is x, 2nd is P(x). We map to y=x/T and I=-(1/T)log P.
# Guard against P<=0: skip those points (1/0 = NaN in gnuplot).
Y(colx) = (colx)/T
IT(colp) = (colp>0 ? (-1.0/T)*log(colp) : 1/0)

# ---- plot ----
plot "gaussian_simulation1.txt" using (Y($1)):(IT($2)) with lines lw 2 title "Gaussian (empirical)", \
     "exponential_shifted.txt" using (Y($1)):(IT($2)) with lines lw 2 title "Exponential-shifted (empirical)", \
     Igauss(x) with lines dt 2 lw 2 title "Gaussian theory: y^2/(4D)"
