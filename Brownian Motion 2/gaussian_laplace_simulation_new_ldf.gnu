D = 1.0
dt = 0.001
Tmax = 200
T = Tmax*dt

set term qt
set xlabel "y = x / T"
set ylabel "I_T(y) = -(1/T) log P_T(x=y T)"
set key left
set grid
Igauss(y) = y*y/(4.0*D)+2.3
Y(colx) = (colx)/T
IT(colp) = (colp>0 ? (-1.0/T)*log(colp) : 1/0)

plot "gaussian_simulation1.txt"  u (Y($1)):(IT($2)) w l lw 2 t "Gaussian (emp)", \
     "exponential_shifted.txt"   u (Y($1)):(IT($2)) w l lw 2 t "Exp-shifted (emp)", \
     "laplace.txt"               u (Y($1)):(IT($2)) w l lw 2 t "Laplace (emp)", \
     Igauss(x) w l dt 2 lw 2 t "Gaussian theory"
