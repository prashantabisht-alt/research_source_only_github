set xlabel "x"
set ylabel "P(x)"
set grid

# Uncomment next line to see tails in log scale
# set logscale y

# Theory parameters
m = 1.0
kBT = 1.0
sigma2 = 2.0 * 1.0 * 10000 * 0.001  # 2 D t  (D=1, t = Tmax*dt)

norm_gauss(x) = 1.0/sqrt(2.0*pi*sigma2) * exp(-x**2/(2.0*sigma2))

# Plot
plot \
  "gaussian_simulation1.txt" using 1:2 with points pt 7 ps 0.3 lc rgb "blue" title "Gaussian Sim (noise ~ N(0,1))", \
  "exponential_simulation1.txt" using 1:2 with points pt 7 ps 0.3 lc rgb "red" title "Exponential Sim (Laplace noise)", \
  norm_gauss(x) with lines lw 2 lc rgb "black" title "Theory Gaussian"
