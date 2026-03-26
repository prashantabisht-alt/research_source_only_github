set logscale xy
set xlabel "t/τ_p"
set ylabel "MSD / (γ^2 v0^2 τ_p^6 / λ^2)"
plot \
  "msd_persistent_norm.txt"    u 1:2 w p lw 2 t "simulation (OU)", \
  "theory_persistent_norm.txt" u 1:2 w p dt 2 lw 2 t "theory (Eq.23)", \
  (1.0/36.0)*x**6 w l dt 3 lw 3 t "t^6 guide", \
  0.1*x**5        w l dt 3 lw 3 t "t^5 guide"
