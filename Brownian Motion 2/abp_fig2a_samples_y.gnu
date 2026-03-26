reset
set term qt
set xlabel "y / σ_y"; set ylabel "PDF"
normal(x)=1.0/sqrt(2.0*pi)*exp(-0.5*x**2)
plot "y_scaled_t0.25.txt" u 1:2 w lp t "t=0.25", \
     "y_scaled_t0.50.txt" u 1:2 w lp t "t=0.50", \
     "y_scaled_t1.00.txt" u 1:2 w lp t "t=1.00", \
     normal(x) w l lw 2 t "N(0,1)"

