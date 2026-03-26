set term qt
set xlabel "z = (v0^2 t^2 - r^2)/(2 DR v0^2 t^3)"
set ylabel "f_r(z)"
plot "fr_t0.25.txt" u 1:2 w lp t "t=0.25", \
     "fr_t0.50.txt" u 1:2 w lp t "t=0.50", \
     "fr_t1.00.txt" u 1:2 w lp t "t=1.00", \
     "fr_theory.txt" u 1:2 w l lw 2 lc rgb "red" t "theory f_r(z)"
