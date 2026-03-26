# Params
D=1.0; dt=0.001; Tmax=200; T=Tmax*dt; pi=3.141592653589793
a = sqrt(2.0*D*dt)

# Empirical helpers
Y(colx)  = (colx)/T
IT(colp) = (colp>0 ? (-1.0/T)*log(colp) : 1/0)

# Gaussian theory (exact finite-T)
Cgauss(T,D) = (1.0/(2.0*T))*log(4.0*pi*D*T)
Igaussian(y)   = y*y/(4.0*D) + Cgauss(T,D)

# Exponential increments: asymptotic + finite-T correction
Iexp_asym(y) = (y>-a/dt ? (1.0/dt)*((dt/a)*y - log(1.0+(dt/a)*y)) : 1/0)
Iexponential(y) = (y>-a/dt ? Iexp_asym(y) + (1.0/(2.0*T))*log( 2.0*pi*T * ((2.0*D + a*y)**2 / (2.0*D)) ) : 1/0)

# Plot
set term qt size 1100,650
set title sprintf("Scaled rate function I_T(y) at T=%.3g", T)
set xlabel "y = x / T"; set ylabel "I_T(y) = -(1/T) ln P_T(x=yT)"
set grid back lw 1; set border lw 1.8; set tics out; set key left box opaque width -2 spacing 1.2
set samples 800; set xrange [-13:13]; set yrange [0:45]

set style line 1 lc rgb "#1f77b4" lw 2.2    # Gaussian (emp)
set style line 2 lc rgb "#d62728" lw 2.2    # Exp (emp)
set style line 4 lc rgb "#000000" lw 2.8 dt 2
set style line 6 lc rgb "#228833" lw 2.8 dt 2

plot \
  "gaussian_simulation1.txt"  u (Y($1)):(IT($2)) w l ls 1 t "Gaussian (emp)", \
  "exponential_shifted.txt"   u (Y($1)):(IT($2)) w l ls 2 t "Exp-shifted (emp)", \
  Igaussian(x)                   w l ls 4 t "Gaussian theory", \
  Iexponential(x)                w l ls 6 t "Exponential theory"
