set xlabel "a1 = (v0 t - x)/(v0 D_R t^2)"
set ylabel "(v0 D_R t^2) P(x,t)  (collapsed density)"
set key top right

plot \
  "a1_edge_t0.25.txt" u 1:2 w l lw 2 t "t=0.25", \
  "a1_edge_t0.50.txt" u 1:2 w l lw 2 t "t=0.50", \
  "a1_edge_t1.00.txt" u 1:2 w l lw 2 t "t=1.00"
