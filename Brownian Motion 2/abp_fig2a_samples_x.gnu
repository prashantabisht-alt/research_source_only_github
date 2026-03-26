reset
set term qt
set xlabel "a1 = (v0 t - x)/(v0 DR t^2)"; set ylabel "PDF"
plot "a1_edge_t0.25.txt" u 1:2 w lp t "t=0.25", \
     "a1_edge_t0.50.txt" u 1:2 w lp t "t=0.50", \
     "a1_edge_t1.00.txt" u 1:2 w lp t "t=1.00", \
     "fx_theory.txt"     u 1:2 w l  lw 2 t "theory f_x"
