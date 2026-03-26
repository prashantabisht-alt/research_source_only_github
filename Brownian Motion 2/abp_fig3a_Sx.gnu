set logscale x
set logscale y
set xlabel "t"
set ylabel "S_x(t; x_0)"
set key bottom left

# ----- simulation parameters -----
pi = 3.141592653589793
v0 = 1.0
x0 = 0.1

# theoretical asymptotic prefactor
Sx_theory(t,DR) = (x0*sqrt(2.0*DR)/(sqrt(pi)*v0)) * t**(-0.5)

# arbitrary slope guide (paper style), anchored at a point (t0,S0)
t0 = 10.0
S0 = 0.1
SlopeGuide(x) = S0 * (x/t0)**(-0.5)*10
# ---- plot your data ----
plot \
 "Sx_DR0.0050.txt" u 1:2 w l lw 2 t "D_R = 0.005", \
 "Sx_DR0.0100.txt" u 1:2 w l lw 2 t "D_R = 0.01",  \
 "Sx_DR0.0500.txt" u 1:2 w l lw 2 t "D_R = 0.05",  \
 Sx_theory(x,0.005) w l dt 2 lw 2 t "theory prefactor (0.005)", \
 Sx_theory(x,0.010) w l dt 2 lw 2 t "theory prefactor (0.01)",  \
 Sx_theory(x,0.050) w l dt 2 lw 2 t "theory prefactor (0.05)",  \
 SlopeGuide(x)      w l dt 3 lw 2 t "paper-style slope ~ t^{-1/2}"
