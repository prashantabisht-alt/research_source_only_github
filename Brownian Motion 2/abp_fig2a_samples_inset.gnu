set xlabel "y / σ_y"
set ylabel "σ_y P(y,t)  (collapsed density)"
set key top right

plot \
  "y_scaled_t0.25.txt" u 1:2 w l lw 2 t "t=0.25", \
  "y_scaled_t0.50.txt" u 1:2 w l lw 2 t "t=0.50", \
  "y_scaled_t1.00.txt" u 1:2 w l lw 2 t "t=1.00", \
  1.0/sqrt(2.0*pi)*exp(-x**2/2.0) w l dt 2 lw 2 t "N(0,1)"
