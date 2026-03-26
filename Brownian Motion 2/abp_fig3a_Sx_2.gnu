set logscale x
set logscale y
set xlabel "u = t D_R"
set ylabel "S_x(u)"
set key bottom left

pi = 3.141592653589793
v0 = 1.0
x0 = 0.1
Sx_u(u) = (x0*sqrt(2.0)/(sqrt(pi)*v0)) * u**(-0.5)

plot \
 "Sx_collapsed_DR0.0050.txt" u 1:2 w l lw 2 t "D_R = 0.005", \
 "Sx_collapsed_DR0.0100.txt" u 1:2 w l lw 2 t "D_R = 0.01",  \
 "Sx_collapsed_DR0.0500.txt" u 1:2 w l lw 2 t "D_R = 0.05",  \
 Sx_u(x) w l dt 3 lw 1.8 t sprintf("const · u^{-1/2} (%.3g)", x0*sqrt(2)/(sqrt(pi)*v0))
